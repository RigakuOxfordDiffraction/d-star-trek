#ifndef DT_CREFLNLIST_H
#define DT_CREFLNLIST_H
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
// Creflnlist.h    Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Creflnlist
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
#include "dtreksys.h"
#include "Cstring.h"
#include "Crefln.h"
#include "Ccrystal.h"  // Needed to reduce reflnlist to asymmetric unit
#include "ResoBins.h"

//+Definitions and constants
// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eReflnFieldType {
  eReflnField_unknown_type,
  eReflnField_int_type,
  eReflnField_float_type,
  eReflnField_Cstring_type
};

enum eReflnlist_states {
  eReflnlist_unknown_state,
  eReflnlist_available_state,
  eReflnlist_pickedpeaks_state,
  eReflnlist_predicted_state,
  eReflnlist_integrating_state,
  eReflnlist_integrated_state,
  eReflnlist_scaled_state,
  eReflnlist_merged_state,
  eReflnlist_averaged_state
};
enum eReflnoutput {
    eReflnoutput_binary,
    eReflnoutput_text,
    eReflnoutput_default
};
enum eLogoutput {
    eLogoutput_log,
    eLogoutput_nolog,
    eLogoutput_default
};

// For table output.

enum { 
        eTableResoStartEnd =                        0x0001,
        eTableIoverSigRange =                       0x0002,
        eTableNumReflns =                           0x0004,
        eTableNumRejects =                          0x0008,
        eTableIoverSig =                            0x0010,
        eTableMaxReso =                             0x0020,
        eTableAvgBackground =                       0x0040,
        eTableAvgBackgroundOverAvgBackgroundSig =   0x0080,
        eTableSharpness =                           0x0100
};

enum {
    eTableVsResolution = 1,
    eTableVsIoverSig = 2
};

//////////////////////////////////////////////////////////////////////////
typedef struct _tagSTRATEGY_REFLN_OVERLAP_CHECK
{
    double      m_dDistTest;  // test 2D center coordinate difference (pix)
    double      m_dOscTest_1; // test oscillation width (deg)
    double      m_dOscStep;   // test oscillation width step (deg)
    double      m_dOscTest_2; // test oscillation width max (deg)

    _tagSTRATEGY_REFLN_OVERLAP_CHECK()
    {
        m_dDistTest  = -1.0;  // default: use the spot size from the header  
        m_dOscTest_1 = 0.5;  
        m_dOscStep   = 0.0;  
        m_dOscTest_2 = m_dOscTest_1; 
    }

    bool bParseInput(Cstring& sArgs)
    {
        std::vector<double>     aParams;
        sArgs.nListToVector(aParams, " ");
        int     nParamsCount = (int)aParams.size();
                                                                    
        switch( nParamsCount )
        {
        case 0:
            break;  // no parameters supplied
        case 1:
            m_dDistTest  = aParams[0];
            break;
        case 2:
            m_dDistTest  = aParams[0];
            m_dOscTest_1 = aParams[1]; 
            m_dOscStep   = 0.0; 
            m_dOscTest_2 = m_dOscTest_1;
            break;
        case 4:
            m_dDistTest  = aParams[0];
            m_dOscTest_1 = aParams[1]; 
            m_dOscStep   = aParams[2]; 
            m_dOscTest_2 = aParams[3];
            break;
        default:
            return false;
        }
	
	    return true;
    }
    
    void vPrint()
    {
        printf("\nStrategy reflection overlap check parameters:\n");
        printf("2D center coordinate difference: %.0f pix\n",   m_dDistTest);
        printf("First rotation width to test: %.2f deg\n",   m_dOscTest_1);
        printf("Last rotation width to test:  %.2f deg\n",   m_dOscTest_2);
        printf("Rotation width step:          %.2f deg\n\n", m_dOscStep);
    }
}STRATEGY_REFLN_OVERLAP_CHECK;
////////////////////////////////////////////////////////////////////////////////////////////

//+Code begin

class Crefln;                  // Forward reference of class Crefln
struct tagBinary_Ref_Header;   // Forward reference of binary header structure.

class DTREK_EXPORT Creflnlist 
{

  friend class Crefln;

protected:
  eReflnlist_states
           m_eThe_State;       // The state of the reflection list

  Cstring  m_sChildLog;        // The history log of each reflection list that was read in.  When we write out, all levels
                               // are incremented by one.
  Cstring  m_sCommentLog;      // Any optional comments added are placed here.
  Cstring  m_sReflnlist_key;   // Unique Reflnlist key for database lookup

  Cstring  m_sEmbeddedHeader;  // Embedded header information.
  int      m_nNumReflns;       // Number of reflections in list
  int      m_nNumReflnsAvail;  // Normally set to 0, and is only calculated in nRead() if set to -1.


  int      m_nNumAlloc;        // Allocated size of m_ppoTheReflns
  int      m_nAllocSize;       // Size to allocate at a single time

  Cstring *m_psIntNames;       // Property names of the Creflns' nInt fields

  Cstring *m_psFloatNames;     // Property names of the Creflns' nFloat fields

  Cstring *m_psCstringNames;   // Property names of the Creflns' nCstring fields

  Crefln  *m_poReflnMin;       // A dummy reflection with minimum values of
                               //    all reflection properties
  Crefln  *m_poReflnMax;       // A dummy reflection with maximum values of
                               //    all reflection properties

  int      m_nNoWriteHeader;   // Do not write header if this var non-zero.
  int     *m_pnFortranFormat;  // Fortran format if specified.
  int     *m_pnFortranFormat2;

  eReflnoutput m_eWriteBinary; // Should we write the output in binary format, text format or the default format?
  eLogoutput m_eWriteLog;      // Should we output log information in the header files?
  tagBinary_Ref_Header* m_poBinaryHeader;
                               // This is the header for a binary reflection list.
                               // Once we have started writting a reflection list, this will be NULL or pointing to an object
                               // Crefln::vWrite() looks at whetherthis is NON-NULL to see if it should write binary or text.
  int   m_nBytesOrLinesReadInHeader;   
                              // Number of bytes or lines read in the header.  This information is used by some routines.


  
  int* m_pnIntFieldsRelative;      // These are loaded when a call to nInsertListFrom() is made.  They allow us to use nInsertReflnFrom()
  int* m_pnFloatFieldsRelative;
  int* m_pnCstringFieldsRelative;


public:
        bool bNeedsSaving(void);
        void NeedsSaving(bool NewVal);
        bool m_bOverwrite;

        Cstring m_csInColor;
        Cstring m_csBeforeColor;
        Cstring m_csAfterColor;
        Cstring m_csCalShape;
        Cstring m_csObsShape;
        bool m_bCentroidDiff;
        bool m_bNeedsSaving;
        int m_nCalculatedSpotsExist;
        int m_nObservedSpotsExist;

  int     *m_pnSortIndex;      // Internal index array for sorts.
  int      m_nIntReflnFields;  // Number of integer fields in *ppoTheReflns...
  int      m_nFloatReflnFields;// Number of float fields in *ppoTheReflns...
  int      m_nCstringReflnFields;// Number of Cstring fields in *ppoTheReflns...

  int      m_nTotalFieldsPlus1;// (Sum of previous 3 variables) + 1

  char    *m_pcSelect;         // Selection field for use in select or exclude

  Crefln **m_ppoTheReflns;     // Pointer to pointers of objects Crefln

  Crefln **m_ppoFreeReflns;    // Pointer to pointers of objects Crefln
  int      m_nNumFree;         // Number used in the free list
  int      m_nNumAllocFree;    // Allocated size of m_ppoFreeReflns
  int      m_nVerbose;         // Verbose flag for level of output

// Some public static Cstrings of common refln field names.
// We have them here as static, so that they need not be duplicated
// in any other class.
// The naming convention for these variables are 'ms_s' for static Cstring,
// then each n, f, s or int, float and string respectively.

DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_snDetNum;
DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_snNonunfFlag;
  static Cstring ms_snH;
  static Cstring ms_snK;
  static Cstring ms_snL;
  static Cstring ms_sfIntensity;
  static Cstring ms_sfSigmaI;
  static Cstring ms_sfProfitIntensity;
  static Cstring ms_sfProfitSigmaI;
  static Cstring ms_snPartialFlag;
  static Cstring ms_sfIntensityPlus;
  static Cstring ms_sfSigmaIPlus;
  static Cstring ms_sfIntensityMinus;
  static Cstring ms_sfSigmaIMinus;
  static Cstring ms_sfObsPx0;
  static Cstring ms_sfObsPx1;
  static Cstring ms_sfObsPxPeak0;
  static Cstring ms_sfObsPxPeak1;
  
  static Cstring ms_sfObsSharpness;

  static Cstring ms_sfObsPeakAreaCount;
  static Cstring ms_sfObsPeakBorderCount;

  static Cstring ms_sfObsRotMid;
  static Cstring ms_sfObsRotEnd;
  static Cstring ms_sfObsRotWidth;
  static Cstring ms_sfObsRotSigma;
  static Cstring ms_sfObsImgMid;
  static Cstring ms_sfObsXmm;
  static Cstring ms_sfObsYmm;
  static Cstring ms_sfObsZmm;
  static Cstring ms_sfCalcPx0;
  static Cstring ms_sfCalcPx1;
  static Cstring ms_sfCalcXmm;
  static Cstring ms_sfCalcYmm;
  static Cstring ms_sfCalcZmm;
  static Cstring ms_sfCalcRotStart;
  static Cstring ms_sfCalcRotMid;
  static Cstring ms_sfCalcRotEnd;
  static Cstring ms_sfCalcRotWidth;
  static Cstring ms_sfCalcMosCoeffA;
  static Cstring ms_sfCalcMosCoeffB;
  static Cstring ms_sfPartiality;
  static Cstring ms_sfResolution;
  static Cstring ms_sfDetResolution;
  static Cstring ms_sfRecipCoord0;
  static Cstring ms_sfRecipCoord1;
  static Cstring ms_sfRecipCoord2;
  static Cstring ms_sfRecipCoordD0;
  static Cstring ms_sfRecipCoordD1;
  static Cstring ms_sfRecipCoordD2;
  static Cstring ms_snReflnNum1;
  static Cstring ms_snReflnNum2;
  static Cstring ms_snOrignlReflnNum;
  static Cstring ms_snRefineBatch;
  static Cstring ms_snDiffFreq;
  static Cstring ms_sfFloatH;
  static Cstring ms_sfFloatK;
  static Cstring ms_sfFloatL;
  static Cstring ms_sfHKLResid;
  static Cstring ms_snTwinID;
  static Cstring ms_sfGonio1;
  static Cstring ms_sfGonio2;
  static Cstring ms_sfGonio3;
  static Cstring ms_sfGonio4;
  static Cstring ms_sfGonio5;
  static Cstring ms_snGonioRotAxis;
  static Cstring ms_sfSvec0;
  static Cstring ms_sfSvec1;
  static Cstring ms_sfSvec2;
  static Cstring ms_sfS0vec0;
  static Cstring ms_sfS0vec1;
  static Cstring ms_sfS0vec2;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidA00;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidA01;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidA11;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidb0;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidb1;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipsoidc;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipseMajorAxis;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipseMinorAxis;
  DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_sfEllipseAxisMajorOffset;
  static Cstring ms_snSourceRefNum;
  static Cstring ms_snCompeteRefNum;


  static Cstring ms_sfPolarz;
  static Cstring ms_sfLorentz;
  static Cstring ms_sfOblique;
  static Cstring ms_sfDeltaPx0;
  static Cstring ms_sfDeltaPx1;
  static Cstring ms_sfDeltaRot;
  static Cstring ms_sfBackground;
  static Cstring ms_sfBackgroundSigma;
  static Cstring ms_ssRejectString;

  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snPackedHKL;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snReducedH;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snReducedK;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snReducedL;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snFplusminusFlag;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_snCentPhase;
  static Cstring ms_ssIsEmpty;

  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_ssBatch;
  DTREK_WIN_DLL_DATA_EXPORT static Cstring ms_ssTag;

  static Cstring ms_sbOvlpFlg;

// Now for some Field indices which go with the above

  int m_nFI_nDetNum;
  int m_nFI_nNonunfFlag;
  int m_nFI_nH;
  int m_nFI_nK;
  int m_nFI_nL;
  int m_nFI_fIntensity;
  int m_nFI_fSigmaI;
  int m_nFI_fProfitIntensity;
  int m_nFI_fProfitSigmaI;
  int m_nFI_nPartialFlag;

  int m_nFI_fIntensityPlus;
  int m_nFI_fSigmaIPlus;
  int m_nFI_fIntensityMinus;
  int m_nFI_fSigmaIMinus;
  int m_nFI_fObsPx0;
  int m_nFI_fObsPx1;
  int m_nFI_fObsPxPeak0;
  int m_nFI_fObsPxPeak1;
  
  int m_nFI_fObsSharpness;
  int m_nFI_fObsPeakAreaCount;
  int m_nFI_fObsPeakBorderCount;

  int m_nFI_fObsPx0Sig;
  int m_nFI_fObsPx1Sig;
  int m_nFI_fObsRotMid;
  int m_nFI_fObsRotEnd;
  int m_nFI_fObsRotWidth;
  int m_nFI_fObsImgMid;
  int m_nFI_fObsXmm;
  int m_nFI_fObsYmm;
  int m_nFI_fObsZmm;
  int m_nFI_fCalcPx0;
  int m_nFI_fCalcPx1;
  int m_nFI_fCalcXmm;
  int m_nFI_fCalcYmm;
  int m_nFI_fCalcZmm;
  int m_nFI_fCalcRotStart;
  int m_nFI_fCalcRotMid;
  int m_nFI_fCalcRotEnd;
  int m_nFI_fCalcRotWidth;
  int m_nFI_fCalcMosCoeffA;
  int m_nFI_fCalcMosCoeffB;
  int m_nFI_fObsRotSigma;
  int m_nFI_fPartiality;
  int m_nFI_fResolution;
  int m_nFI_fDetResolution;
  int m_nFI_fRecipCoord0;
  int m_nFI_fRecipCoord1;
  int m_nFI_fRecipCoord2;
  int m_nFI_fRecipCoordD0;
  int m_nFI_fRecipCoordD1;
  int m_nFI_fRecipCoordD2;
  int m_nFI_nReflnNum1;
  int m_nFI_nReflnNum2;
  int m_nFI_nOrignlReflnNum;
  int m_nFI_nRefineBatch;
  int m_nFI_nDiffFreq;
  int m_nFI_fFloatH;
  int m_nFI_fFloatK;
  int m_nFI_fFloatL;
  int m_nFI_fHKLResid;
  int m_nFI_fPolarz;
  int m_nFI_fLorentz;
  int m_nFI_fOblique;
  int m_nFI_fDeltaPx0;
  int m_nFI_fDeltaPx1;
  int m_nFI_fDeltaRot;
  int m_nFI_fBackground;
  int m_nFI_fBackgroundSigma;
  int m_nFI_nTwinID;
  int m_nFI_sRejectString;
  int m_nFI_fGonio1;
  int m_nFI_fGonio2;
  int m_nFI_fGonio3;
  int m_nFI_fGonio4;
  int m_nFI_fGonio5;
  int m_nFI_nGonioRotAxis;

  int m_nFI_fSvec[3];
  int m_nFI_fS0vec[3];
  int m_nFI_fEllipsoidA00;
  int m_nFI_fEllipsoidA01;
  int m_nFI_fEllipsoidA11;
  int m_nFI_fEllipsoidb0;
  int m_nFI_fEllipsoidb1;
  int m_nFI_fEllipsoidc;
  
  int m_nFI_fEllipseMajorAxis;
  int m_nFI_fEllipseMinorAxis;
  
  int m_nFI_fEllipseAxisMajorOffset;
  
  int m_nFI_nSourceRefNum;
  int m_nFI_nCompeteRefNum;
  int m_nFI_nPackedHKL;
  int m_nFI_nReducedH;
  int m_nFI_nReducedK;
  int m_nFI_nReducedL;
  int m_nFI_nFplusminusFlag;
  int m_nFI_nCentPhase;
  int m_nFI_sBatch;
  int m_nFI_sTag;

  int m_nFI_bOvlpFlg;

// Some static strings for the reflection selection mechanism

  static Cstring ms_sOpEqual;
  static Cstring ms_sOpLessEqual;
  static Cstring ms_sOpGreaterEqual;
  static Cstring ms_sOpNotEqual;
  static Cstring ms_sOpLess;
  static Cstring ms_sOpGreater;
  static Cstring ms_sOpAnd;
  static Cstring ms_sOpOr;

  static Cstring ms_sOpEqual_Fortran;
  static Cstring ms_sOpLessEqual_Fortran;
  static Cstring ms_sOpGreaterEqual_Fortran;
  static Cstring ms_sOpNotEqual_Fortran;
  static Cstring ms_sOpLess_Fortran;
  static Cstring ms_sOpGreater_Fortran;
  static Cstring ms_sOpAnd_Fortran;
  static Cstring ms_sOpOr_Fortran;

  static Cstring ms_sOpASMD;
  static Cstring ms_sOpAddEqual;
  static Cstring ms_sOpSubEqual;
  static Cstring ms_sOpMulEqual;
  static Cstring ms_sOpDivEqual;
  static Cstring ms_sOpSetEqual;
  
  static Cstring ms_asTableHeaders[];
  static Cstring ms_asTableFormat[];

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

////////////////////////////////////////////////////////////////////////

  Creflnlist ();                     // Construct an empty Reflnlist
  Creflnlist (const Cstring& sFilename, const int nMaxToRead=-1, const int nFirstToRead = 0, const bool bIgnoreReadError = FALSE); // Construct a Reflnlist from a file
  Creflnlist (Creflnlist& oReflnlist, bool bSelect = FALSE);
                                     // Clone a Reflnlist, possibly from
                                     //   selected Reflns in another list
  Creflnlist (const int nNum);       // Construct an empty Reflnlist with space
                                     //   for nNum reflections

  Creflnlist (const int nIntFields, const Cstring* psIntNewNames,
              const int nFloatFields, const Cstring* psFloatNewNames,
              const int nCstringNewFields,
              const Cstring* psCstringNewNames,
              const int nNum);


  ~Creflnlist ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes (none are virtual and cannot be overridden!?!
////////////////////////////////////////////////////////////////////////

void LoadViewOptions(Cstring csInColor, Cstring csBeforeColor,
                        Cstring csAfterColor, Cstring csObsShape, Cstring csCalShape);
Crefln* poGetClosestRefln(float x, float y,int nIndexX = -1,int nIndexY = -1);
float m_fPercentRefToDisplay;

void vAddSpot(float x, float y, float intensity, float sigma, float reso,
                         float RotMid, float RotWidth);
int nDeleteSpot(float fPx0, float fPx1);

int  nInitValues(void);
inline int nGetNumReflns(void) { return (m_nNumReflns); }
inline int nGetNumFields(void) { return (m_nTotalFieldsPlus1 - 1); }

inline Crefln* poGetRefln(const int nRefNum) 
                              { return (m_ppoTheReflns[nRefNum]); }

inline Crefln* poGetReflnSafe(const int nRefNum) 
{ 
    if( nRefNum >= m_nNumReflns || nRefNum < 0 )
        return NULL;
    
    return (m_ppoTheReflns[nRefNum]);
}


Crefln& operator[] (const int nIndex) { return *poGetRefln(nIndex); };
int  nSelectField(const eReflnFieldType eType,
                  const int nSelectField, const char cValue='\0',
                  char *pcSelectIn = NULL);

int  nSelect(const int nFI, const Cstring& sOperation,
             const float fTestValue, const bool bIncExc=FALSE);
int  nSelect(const int nFI, const Cstring& sOperation,
             const int nTestValue, const bool bIncExc=FALSE);
int  nSelect(const int nFI, const Cstring& sOperation,
             const Cstring& sTestValue, const bool bIncExc=FALSE);

int  nSelect(const Cstring& sCompareString);

int  nInsert(Crefln* poRefln, int nIndex = -1);
int  nInsert(const int nH, const int nK, const int nL);
int  nInsertSorted(Crefln* poRefln,const eReflnFieldType eType,const int nField);

int  nDelete(const int nReflnNum);
int  nDelete(const int nFlag, int *pnArray);
int  nDelete(const char *pcSelectIn=NULL, const bool bIncExc=FALSE);
int  nDeleteSorted(const int nReflnNumSorted,int* pnArray1 = NULL,int* pnArray2 = NULL,int* pnArray3 = NULL);
int  nDeleteBatchID(Cstring& rsBatchIDTemplate);
void vDeleteAll(void);

#define DTREK_REFLNLIST_DL_RANK          0x0001
#define DTREK_REFLNLIST_DL_MARK_ONLY     0x0002
#define DTREK_REFLNLIST_DL_FROM_FIND     0x0004
#define DTREK_REFLNLIST_DL_FROM_REFINE   0x0008
int  nDeleteRing(double* afIceRingReso, int nRings, double fFull2ThetaWidthToReject, double fWavelength, DTREK_WORD wCtrl=0U);

eReflnoutput eSetWriteBinary(eReflnoutput eSet) { eReflnoutput ePrev = m_eWriteBinary; m_eWriteBinary = eSet; return ePrev;};
int  nSetWriteLog(eLogoutput eSet) { m_eWriteLog = eSet; return 0;};
int  nRead(const Cstring& sFilename, int nMaxToRead=-1, int nFirstToRead = 0, bool bIgnoreReadError = FALSE);
static int  nGetNumReflns(const Cstring& sFilename);

inline void vSetOverwrite(const bool bFlag) { m_bOverwrite = bFlag; }
int  nWrite(const Cstring& sFilename, int *pnIndex = NULL, 
            const char *pcSelectIn = NULL, FILE* pFOutUser = NULL,
            const int nFirst=-1, const int nLast=-1);
void vWritePartial(const Cstring& rsFilename);
int nUpdateCrystal(Cstring& sFileName,Ccrystal* poCrystal = NULL,Cstring* psCommentLog = NULL);

void vSort(const eReflnFieldType eType, const int nField, int* pnIndex);
void vSort2(const eReflnFieldType eType1, const int nField1,const eReflnFieldType eType2, const int nField2,int* pnIndexIn);

int  nFindFirst(int nField,const int nValue,int* pnIndex = NULL,bool bFindFirstGreater = FALSE);                 // Finds the first entry (relative to pnIndex) that agrees with the specified INTEGER value.
int  nFindFirst(int nField,const double fValue,double fTol,int* pnIndex = NULL,bool bFindFirstGreater = FALSE);  // Finds the first entry (relative to pnIndex) that agrees within tolerance with the specified FLOATING POINT value.
int  nFindFirst(int nField,const Cstring& sValue,int* pnIndex = NULL,bool bFindFirstGreater = FALSE);            // Finds the first entry (relative to pnIndex) that matches the specified STRING.


int nList(const int nNum = 0);
int nListFields(void);
int nList(Cimage_header& oHeader);

int nGetFieldIndex(const Cstring& sFieldName,
                   const eReflnFieldType eType = eReflnField_unknown_type);
Cstring sGetFieldName(const int nFieldIndex, const eReflnFieldType eType) const;

int nExpandGetField(const Cstring&  sFieldName,
                    const eReflnFieldType eType = eReflnField_unknown_type);

int nExpandRefln(const Creflnlist& oOther, const char *pcSelectIn=NULL);

int nExpandRefln(const int nIntNewFields, // Should have to give names of fields
                 const Cstring* psIntNewNames,
                 const int nFloatNewFields = 0,
                 const Cstring* psFloatNewNames = NULL,
                 const int nCstringNewFields = 0,
                 const Cstring* psCstringNewNames = NULL);

int nInsertListFrom(Creflnlist& oReflnlist, const char *pcSelectIn=NULL, const int * pnSortOrder=NULL,int nNumToInsert = -1);
int nInsertReflnFrom(Crefln& oReflnNative,Crefln& oReflnOther);


int nConvertFromAnom(void);
int nCalcRecipLatticePoints(Cimage_header& oHeader,double& fMinTwinIDTolerance);

#define DTREK_REFLNLIST_RT_RANK          0x0001
int nResoTable(double         fResoMin = -1.0, 
               double         fResoMax = -1.0,
               int            nNumberOfResoBins = 10,
               int            nNumberOfIOverSigmaBins = 10,
               int            nTableType = eTableVsResolution,
               int            nColsToPrint = -1,
               Cimage_header* poHeader = NULL,
               DTREK_WORD     wCtrl = 0U
                );

void vRank(Cimage_header* poHeader, 
           double fResoMin,      
           double fResoMax,      
           int nNumberOfResoBins,
           int nMinReflnsForStats,
           bool bEvalExpTime=false);

bool bEvaluateExpTime(Cimage_header* poHeader, double fResoMax, double fTargetIOverSigma, double fIOverSigmaTimePower);


int nReduce(Ccrystal &roCrystal, const int nAnomFlag=0,const int nCentricCheck=1) 
            { return nReduce(*roCrystal.m_poSpacegroup, nAnomFlag, nCentricCheck); };
int nReduce(Cspacegroup &roSpace, const int nAnomFlag=0,
                                 const int nCentricCheck=1);
int nReindex(const float *pfMatIn);
int nReindex(const double *pfMatIn);
int nAddResol(Ccrystal &oCrystal);
int nAddDetResol(Cimage_header *poHeader);


inline int *pnGetSortIndex(void) { return (m_pnSortIndex); }
inline  void vSetFortranFormat(int* pnFortranFormat,int* pnFortranFormat2) { m_pnFortranFormat = pnFortranFormat;m_pnFortranFormat2 = pnFortranFormat2; };
inline  void vSetNoWriteHeader(void)
     { m_nNoWriteHeader = 1; }
inline  void vSetWriteHeader(void)
     { m_nNoWriteHeader = 0; }
inline bool bIsAvailable(void)
  {  return (m_eThe_State != eReflnlist_unknown_state); }
inline void vSetState(const eReflnlist_states eState) { m_eThe_State = eState; }

int nOverlapCheck(const int nFI_nHead, const int nFI_nNext, double fRotWidth, double fScanRotStart=-999.0);
int nOverlapCheck(double dSpotIntegrBoxWidth, double fRotWidth, double  fScanRotStart, bool bDisplay=false);
bool bTwoReflnsOnSameFrame(int index1, int index2, double fScanRotStart, double fRotWidth);

void vSortAbsFloatHeap(const int nField, int* pnIndex);
int  nExpand(const int nAdditional);
int nDiffs(int nNumConstantFields,eReflnFieldType* peConstantFields,int* pnConstantFields,float* pfConstantFieldTolerances);

int nDisplayLog(int nVerbose);
int nLogComment(const Cstring& sComment,bool bAdd);
Cstring sGetComment() { return m_sCommentLog; };


// Embedded functions.
Ccrystal* poGetCrystal();
Cimage_header* poGetHeader();
int nPutHeaderInfo(Cimage_header& oHeader,Cstring sPattern,bool bDeleteOld);
int nPutCrystal(Ccrystal& oCrystal);
int nPutCrystal(Cimage_header& oHeader) { Ccrystal oCrystal(oHeader); if (oCrystal.bIsAvailable()) return nPutCrystal(oCrystal); else return 1; };
int nAddExtra(Creflnlist& oList);


protected:

    int nWriteHeader(FILE* pFOut, const char *pcSelectIn = NULL,const bool bBinary = FALSE);
    int nReadHeader(const Cstring& sFilename, char** ppcBuffer, int nMaxCharsInBuffer, char** pcExtraData, tagBinary_Ref_Header& oHeader, int& nBytesOrLinesRead);

    // Logging functions.
    int nEraseExtra();
    int nParseExtra(char* pcInputLog);
    int nGenerateLog(Cstring& sLogbuf,const Cstring& sOutputname,tagBinary_Ref_Header&);


    void vSetFlag(const int nIndex, const int nResultField, const int nResultMask);

    void vSortFloatHeap(const int nField, int* pnIndex,int nFirstRefln=0,int nNumReflns=0,bool bUseIndexOrder=FALSE);
    void vSortIntHeap(const int nField, int* pnIndex,int nFirstRefln=0,int nNumReflns=0,bool bUseIndexOrder=FALSE);
    void vSortStringHeap(const int nField, int *pnIndex,int nFirstRefln=0,int nNumReflns=0,bool bUseIndexOrder=FALSE);
    void vSortIntRadix(const int nField, int *pnIndex,int nFirstRefln=0, int nNumReflns=0);
    void vSortStringRadix(const int nField, int *pnIndex,int nFirstRefln=0, int nNumReflns=0);
    void vSortStringRadixRightAdjusted(const int nField, int *pnIndex, int nFirstRefln=0, int nNumReflns=0, bool bUseInputIndexArray=false); 

    virtual void vUpdateFieldIndices(void);
    int  nDeleteFree(void);

public:
    void vAnalyzeBkgVsResolution();
    void vObliqueIncidenceCorrection(const int nApply = 0, Cimage_header *poHeader=NULL);

    void vSetExternalSortIndex(int* pnSortIndex);

    void vRemoveReflnsInResolutionRingAreas(Cimage_header* poHeader, double dRingWidth, DTREK_WORD wCtrl = 0U);
    void vRemoveOverlapReflns(double fRotWidth);
    void vRemoveIntensityOrSigmaMarkedReflections();


    void vSetColorsForDisplay();

    void vGetResoLimits(double& dResoMin, double& dResoMax, int nResoFieldIndex);
    void vToStratResoBins(CResoStrategyBins* poResoBins);

private:
    typedef enum
    {
        enColorBlue    = 0,
        enColorRed     = 1,
        enColorGreen   = 2,
        enColorMagenta = 3,
        enColorSkip    = 4,
        enColorYellow  = 5
    }enumDisplayColor;

};// end of class Creflnlist


// This type of float is used in the binary output.

typedef float CReffloat;

// Binary header information.

struct tagBinary_Ref_Header {

    char cSizeFloat;    // These two fields need to be chars since they are used BEFORE reversing the order of bytes.
    char cSizeInt;
    int nVersion;       // Version.Probably used if we have more "extra" data, and want to support the extra fields.
    int nNumRefs;
    int nMaxString;     // Maximum size of a string (including NULL terminator).  Placed here for core dump saftey protection.

    int nExtraSize;     // Size of any extra information following Binary_Ref_Header that must be read in.

    int a3nFields[3];

    union {
        float _fTemp;
        struct  {
            char a16cSampleFloat[16];
            char a16cSampleInt[16];
        } _structbrace ;
    } _unionbrace ;

    int nLoadHeader(Creflnlist& oReflnlist,const char* pcSelectIn=NULL);    // only nExtraSize is not filled in.
    int nWriteHeader(FILE* pFOut);                                          // After the nExtraSize field is filled in, then this can be called.
    int nReadHeader(const Cstring& sFilename, char** ppcBuffer, char** ppcExtra, int nMaxSearch, int& nCharsRead);

    Crefln& oLoadRefln(Crefln& oRefln,char* pcRecord);  // Changes endian-ness in fields so that they are consistent with record fields.
                                                        // Then loads the given reflection.
    int nRead(char* pcBuf,Crefln& oRefln);
    int nWrite(FILE* pFOut,Crefln& oRefln,const char* pcSelectIn=NULL);

    static char* ms_pcId;
    static CReffloat ms_fCodeFloat;     // Code float written into the data file.
    static int ms_nCodeInt;             // Code integer written into the data file.
    static bool  ms_bVax2IEEE;          // Set when nReadHeader() called to indicate the types of conversions necc.
    static bool  ms_bSwapFloat;
    static bool  ms_bSwapLong;
    static Cstring ms_sFields;          // Fields saved in a string.  This is used when writting the log information.

};

#endif   // DT_CREFLNLIST_H
