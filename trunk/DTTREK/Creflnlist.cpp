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
// Creflnlist.cc        Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Creflnlist which implements
//    the reflection list encapsulation of d*TREK.
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

#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include "Dtrek.h"
#include "Creflnlist.h"         // Class definition and prototypes
#include "dtreklogic.h"         // for selection mech comparisons

#include "Cdetector.h"          // Needed for nAddDetResol
#include "Csource.h"            // Needed for nAddDetResol
#include "Crotation.h"
#include "Cgoniometer.h"

#include "ResoBins.h"
#include "dtsvd.h"

#ifdef WIN32
#include <time.h>
#include <io.h>
#endif

//mrp#endif

#include "raxis.h"

using namespace std;

//+Definitions, constants, and initialization of static member variables

#define INITIAL_INT        0
#define INITIAL_FLOAT   -999.0
#define INITIAL_STRING    ""
#define ALLOC_SIZE      10000

Cstring Creflnlist::ms_snDetNum         = "nDetector_number";
Cstring Creflnlist::ms_snNonunfFlag     = "nNonunf_flag";
Cstring Creflnlist::ms_snH              = "nH";
Cstring Creflnlist::ms_snK              = "nK";
Cstring Creflnlist::ms_snL              = "nL";
Cstring Creflnlist::ms_sfIntensity      = "fIntensity";
Cstring Creflnlist::ms_sfSigmaI         = "fSigmaI";
Cstring Creflnlist::ms_sfProfitIntensity= "fIntensityProfit";
Cstring Creflnlist::ms_sfProfitSigmaI   = "fSigmaIProfit";
Cstring Creflnlist::ms_snPartialFlag    = "nPartialFlag";
Cstring Creflnlist::ms_sfIntensityPlus  = "fIntensity+";
Cstring Creflnlist::ms_sfSigmaIPlus     = "fSigmaI+";
Cstring Creflnlist::ms_sfIntensityMinus = "fIntensity-";
Cstring Creflnlist::ms_sfSigmaIMinus    = "fSigmaI-";
Cstring Creflnlist::ms_sfObsPx0         = "fObs_pixel0";
Cstring Creflnlist::ms_sfObsPx1         = "fObs_pixel1";
Cstring Creflnlist::ms_sfObsPxPeak0     = "fObs_pixel_peak0";
Cstring Creflnlist::ms_sfObsPxPeak1     = "fObs_pixel_peak1";

Cstring Creflnlist::ms_sfObsSharpness   = "fObs_sharpness";

Cstring Creflnlist::ms_sfObsPeakAreaCount   = "fObs_peak_area";
Cstring Creflnlist::ms_sfObsPeakBorderCount = "fObs_peak_border_count";

Cstring Creflnlist::ms_sfObsRotMid      = "fObs_rot_mid";
Cstring Creflnlist::ms_sfObsRotEnd      = "fObs_rot_end";
Cstring Creflnlist::ms_sfObsRotWidth    = "fObs_rot_width";
Cstring Creflnlist::ms_sfObsRotSigma    = "fObs_rot_sigma";
Cstring Creflnlist::ms_sfObsImgMid      = "fObs_img_mid";
Cstring Creflnlist::ms_sfObsXmm         = "fObs_Xmm";
Cstring Creflnlist::ms_sfObsYmm         = "fObs_Ymm";
Cstring Creflnlist::ms_sfObsZmm         = "fObs_Zmm";
Cstring Creflnlist::ms_sfCalcPx0        = "fCalc_pixel0";
Cstring Creflnlist::ms_sfCalcPx1        = "fCalc_pixel1";
Cstring Creflnlist::ms_sfCalcXmm        = "fCalc_Xmm";
Cstring Creflnlist::ms_sfCalcYmm        = "fCalc_Ymm";
Cstring Creflnlist::ms_sfCalcZmm        = "fCalc_Zmm";
Cstring Creflnlist::ms_sfCalcRotStart   = "fCalc_rot_start";
Cstring Creflnlist::ms_sfCalcRotMid     = "fCalc_rot_mid";
Cstring Creflnlist::ms_sfCalcRotEnd     = "fCalc_rot_end";
Cstring Creflnlist::ms_sfCalcRotWidth   = "fCalc_rot_width";
Cstring Creflnlist::ms_sfCalcMosCoeffA  = "fCalc_mos_CoeffA";
Cstring Creflnlist::ms_sfCalcMosCoeffB  = "fCalc_mos_CoeffB";
Cstring Creflnlist::ms_sfPartiality     = "fCalc_partiality";
Cstring Creflnlist::ms_sfResolution     = "fResolution";
Cstring Creflnlist::ms_sfDetResolution  = "fDetResolution";
Cstring Creflnlist::ms_sfRecipCoord0    = "fRecip0";
Cstring Creflnlist::ms_sfRecipCoord1    = "fRecip1";
Cstring Creflnlist::ms_sfRecipCoord2    = "fRecip2";
Cstring Creflnlist::ms_sfRecipCoordD0   = "fRecipDeriv0";
Cstring Creflnlist::ms_sfRecipCoordD1   = "fRecipDeriv1";
Cstring Creflnlist::ms_sfRecipCoordD2   = "fRecipDeriv2";
Cstring Creflnlist::ms_snReflnNum1      = "nReflnNum1";
Cstring Creflnlist::ms_snReflnNum2      = "nReflnNum2";
Cstring Creflnlist::ms_snOrignlReflnNum = "nOriginalReflnNum";
Cstring Creflnlist::ms_snRefineBatch    = "nRefineBatch";
Cstring Creflnlist::ms_snDiffFreq       = "nFrequency";
Cstring Creflnlist::ms_sfFloatH         = "fFloatH";
Cstring Creflnlist::ms_sfFloatK         = "fFloatK";
Cstring Creflnlist::ms_sfFloatL         = "fFloatL";
Cstring Creflnlist::ms_sfHKLResid       = "fHKLResid";
Cstring Creflnlist::ms_snTwinID         = "nTwinID";
Cstring Creflnlist::ms_sfGonio1         = "fGonioAxis1";
Cstring Creflnlist::ms_sfGonio2         = "fGonioAxis2";
Cstring Creflnlist::ms_sfGonio3         = "fGonioAxis3";
Cstring Creflnlist::ms_sfGonio4         = "fGonioAxis4";
Cstring Creflnlist::ms_sfGonio5         = "fGonioAxis5";
Cstring Creflnlist::ms_snGonioRotAxis   = "nGonioRotAxis";
Cstring Creflnlist::ms_sfPolarz         = "fCalc_polarz";
Cstring Creflnlist::ms_sfLorentz        = "fCalc_lorentz";
Cstring Creflnlist::ms_sfOblique        = "fCalc_oblique";
Cstring Creflnlist::ms_sfDeltaPx0       = "fDeltaPx0";
Cstring Creflnlist::ms_sfDeltaPx1       = "fDeltaPx1";
Cstring Creflnlist::ms_sfDeltaRot       = "fDeltaRot";
Cstring Creflnlist::ms_sfBackground     = "fBackground";
Cstring Creflnlist::ms_sfBackgroundSigma= "fBackgroundSigma";
Cstring Creflnlist::ms_ssRejectString   = "sRejectString";


Cstring Creflnlist::ms_sfSvec0          = D_K_IntegratefSvec0;
Cstring Creflnlist::ms_sfSvec1          = D_K_IntegratefSvec1;
Cstring Creflnlist::ms_sfSvec2          = D_K_IntegratefSvec2;
Cstring Creflnlist::ms_sfS0vec0         = D_K_IntegratefS0vec0;
Cstring Creflnlist::ms_sfS0vec1         = D_K_IntegratefS0vec1;
Cstring Creflnlist::ms_sfS0vec2         = D_K_IntegratefS0vec2;

Cstring Creflnlist::ms_sfEllipsoidA00   = "fEllipsoidA00";
Cstring Creflnlist::ms_sfEllipsoidA01   = "fEllipsoidA01";
Cstring Creflnlist::ms_sfEllipsoidA11   = "fEllipsoidA11";
Cstring Creflnlist::ms_sfEllipsoidb0    = "fEllipsoidb1";
Cstring Creflnlist::ms_sfEllipsoidb1    = "fEllipsoidb0";
Cstring Creflnlist::ms_sfEllipsoidc    = "fEllipsoidc";
Cstring Creflnlist::ms_sfEllipseMajorAxis         = "fEllipseAxisMajor";
Cstring Creflnlist::ms_sfEllipseMinorAxis         = "fEllipseAxisMinor";
Cstring Creflnlist::ms_sfEllipseAxisMajorOffset   = "fEllipseAxisMajorOffset";

Cstring Creflnlist::ms_snSourceRefNum   = "nSourceRef";
Cstring Creflnlist::ms_snCompeteRefNum  = "nCompeteRef";
Cstring Creflnlist::ms_snPackedHKL      = "nPackedHKL";
Cstring Creflnlist::ms_snReducedH       = "nReducedH";
Cstring Creflnlist::ms_snReducedK       = "nReducedK";
Cstring Creflnlist::ms_snReducedL       = "nReducedL";
Cstring Creflnlist::ms_snFplusminusFlag = "nAnomFlag";
Cstring Creflnlist::ms_snCentPhase      = "nCentPhase";
Cstring Creflnlist::ms_ssIsEmpty        = "?";
Cstring Creflnlist::ms_ssBatch          = "sBatch";
Cstring Creflnlist::ms_ssTag            = "sTag";

Cstring Creflnlist::ms_sOpEqual         = "==";
Cstring Creflnlist::ms_sOpLessEqual     = "<=";
Cstring Creflnlist::ms_sOpGreaterEqual  = ">=";
Cstring Creflnlist::ms_sOpNotEqual      = "!=";
Cstring Creflnlist::ms_sOpLess          = "<";
Cstring Creflnlist::ms_sOpGreater       = ">";
Cstring Creflnlist::ms_sOpAnd           = "&";
Cstring Creflnlist::ms_sOpOr            = "|";

Cstring Creflnlist::ms_sOpEqual_Fortran         =".EQ.";
Cstring Creflnlist::ms_sOpLessEqual_Fortran     =".LE.";
Cstring Creflnlist::ms_sOpGreaterEqual_Fortran  =".GE.";
Cstring Creflnlist::ms_sOpNotEqual_Fortran      =".NE.";
Cstring Creflnlist::ms_sOpLess_Fortran          =".LT.";
Cstring Creflnlist::ms_sOpGreater_Fortran       =".GT.";
Cstring Creflnlist::ms_sOpAnd_Fortran           = ".AND.";
Cstring Creflnlist::ms_sOpOr_Fortran            = ".OR.";

Cstring Creflnlist::ms_sOpASMD          = "+-*/";
Cstring Creflnlist::ms_sOpAddEqual      = "+=";
Cstring Creflnlist::ms_sOpSubEqual      = "-=";
Cstring Creflnlist::ms_sOpMulEqual      = "*=";
Cstring Creflnlist::ms_sOpDivEqual      = "/=";
Cstring Creflnlist::ms_sOpSetEqual      = "=";

typedef bool (* y_prbLOGIC_F_F)(const float fA, const float fB);
typedef bool (* y_prbLOGIC_N_N)(const int nA, const int nB);
typedef bool (* y_prbLOGIC_S_S)(const Cstring sA, const Cstring sB);

//const double      cfIOverSigmaMinimumTargetInTargetResoShell = 10.0;
//const int         cnMinimumAcceptableNumberOfReflnsPerShell = 5;

const double      cfMinimumResolutionForLinearApproximation = 1.5; // angstrom
const int         cnMinimumNumberOfPointsForFitting = 10;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Creflnlist::Creflnlist()
{
 (void) nInitValues();
 m_eThe_State = eReflnlist_available_state;
}

Creflnlist::Creflnlist(const Cstring& sFilename, const int nMaxToRead,
		       const int nFirstToRead, const bool bIgnoreReadError)
{
  int nStat;
  (void) nInitValues();
  nStat = nRead(sFilename,nMaxToRead,nFirstToRead);
  if (0 == nStat) 
      m_eThe_State = eReflnlist_available_state;
}

Creflnlist::Creflnlist (const int nIntFields, const Cstring* psIntNewNames,
            const int nFloatFields, const Cstring* psFloatNewNames,
            const int nCstringNewFields = 0,
            const Cstring* psCstringNewNames = NULL,
            const int nNum = 0)
{
  int i;

  nInitValues();

  m_eThe_State         = eReflnlist_unknown_state;
  m_sReflnlist_key     = "unknown";

  m_pnSortIndex        = NULL;
  m_fPercentRefToDisplay = 1.0;
  m_bNeedsSaving       = 0;

  m_nIntReflnFields    = max(3, nIntFields);   // Must have at least 3 int
  m_nFloatReflnFields  = max(2, nFloatFields); //   and 2 floats
  m_nCstringReflnFields = max(0, nCstringNewFields);
  m_nTotalFieldsPlus1  = m_nIntReflnFields + m_nFloatReflnFields
                                           + m_nCstringReflnFields + 1;
  m_psIntNames        = new Cstring [   m_nIntReflnFields ];
  m_psFloatNames      = new Cstring [ m_nFloatReflnFields ];

  if (0 < m_nCstringReflnFields)
    m_psCstringNames = new Cstring [ m_nCstringReflnFields ];
  else
    m_psCstringNames = NULL;

  m_pcSelect          = new char [m_nTotalFieldsPlus1];
  for (i = 0; i < m_nTotalFieldsPlus1; i++) *(m_pcSelect+i) = (char) 0;

  for (i = 0; i < m_nIntReflnFields; i++)
    *(m_psIntNames+i) = *(psIntNewNames+i);
  for (i = 0; i < m_nFloatReflnFields; i++)
    *(m_psFloatNames+i) = *(psFloatNewNames+i);
  for (i = 0; i < m_nCstringReflnFields; i++)
    *(m_psCstringNames+i) = *(psCstringNewNames+i);

  // Check for the required names!

  if (*m_psIntNames != ms_snH)
    cout << "WARNING in Creflnlist::... first integer field is not "
         << ms_snH << "!\n";
  if (*(m_psIntNames+1) != ms_snK)
    cout << "WARNING in Creflnlist::... second integer field is not "
         << ms_snK << "!\n";
  if (*(m_psIntNames+2) != ms_snL)
    cout << "WARNING in Creflnlist::... third integer field is not "
         << ms_snL << "!\n";

  if (*m_psFloatNames != ms_sfIntensity)
    cout << "WARNING in Creflnlist::... first float field is not "
         << ms_sfIntensity << "!\n";
  if (*(m_psFloatNames+1) != ms_sfSigmaI)
    cout << "WARNING in Creflnlist::... second float field is not "
         << ms_sfSigmaI << "!\n";

  m_nNumReflns        = 0;
  m_nAllocSize        = ALLOC_SIZE;
  m_poReflnMin        = NULL;
  m_poReflnMax        = NULL;

  m_poBinaryHeader = NULL;

  // Allocate some space for pointer to reflection pointers

  if (0 < nNum)
    m_nNumAlloc = nNum;
  else
    m_nNumAlloc = ALLOC_SIZE;

  m_ppoTheReflns = new Crefln* [m_nNumAlloc];
  for (i = 0; i < m_nNumAlloc; i++)
    *(m_ppoTheReflns+i) = NULL;

  m_nNoWriteHeader = 0;
  vUpdateFieldIndices();
}

Creflnlist::Creflnlist(const int nNum)
{
  int nStat;
  (void) nInitValues();
  nStat = nExpand(nNum);  // Allocate memory for at least nNum reflections
  if (0 == nStat) 
      m_eThe_State = eReflnlist_available_state;
  vUpdateFieldIndices();
}

Creflnlist::Creflnlist(Creflnlist& oReflnlistIn, bool bSelect)
{
  // Create a reflection list with the same fields as the input reflection list.
  int nStat;
  (void) nInitValues();
  nStat = nExpandRefln(oReflnlistIn.m_nIntReflnFields-3,   // Less H, K, L
                       &oReflnlistIn.m_psIntNames[3],
                       oReflnlistIn.m_nFloatReflnFields-2, // Less Intensity, SigmaI
                       &oReflnlistIn.m_psFloatNames[2],
                       oReflnlistIn.m_nCstringReflnFields,
                       &oReflnlistIn.m_psCstringNames[0]);

  if (0 == nStat) 
      m_eThe_State = eReflnlist_available_state;

  // For now we can't copy ...
  if (bSelect) m_eThe_State = eReflnlist_unknown_state;
  vUpdateFieldIndices();
  nAddExtra(oReflnlistIn);
}



Creflnlist::~Creflnlist()
{
//  cout << "Creflnlist destructor called!\n";
  int i;

  (void) nDeleteFree();

  for (i = 0; i < m_nNumReflns; i++)
    {  // Delete each reflection
      if (NULL != m_ppoTheReflns[i])
        {
          delete m_ppoTheReflns[i];
          m_ppoTheReflns[i] = NULL;
        }
    }
  if (m_ppoTheReflns)
    {
      delete [] m_ppoTheReflns;
      m_ppoTheReflns = NULL;
    }
  m_nNumReflns = 0;
  m_nNumAlloc  = 0;
  if (NULL != m_poReflnMin)
    {
      delete m_poReflnMin;
      m_poReflnMin = NULL;
    }
  if (NULL != m_poReflnMax)
    {
      delete m_poReflnMax;
      m_poReflnMax = NULL;
    }

  delete [] m_psIntNames;      // Is there a memory leak here since the Cstrings
  m_psIntNames = NULL;
  delete [] m_psFloatNames;    //  themselves are not deleted? Ans: Probably.
  m_psFloatNames = NULL;

  if (NULL != m_psCstringNames)
    {
      delete [] m_psCstringNames;
      m_psCstringNames = NULL;
    }
  delete [] m_pcSelect;
  m_pcSelect = NULL;

  if (NULL != m_pnSortIndex)
    {
      delete [] m_pnSortIndex;
      m_pnSortIndex = NULL;
    }
  if (NULL != m_poBinaryHeader) 
      delete m_poBinaryHeader;

  if (NULL != m_pnIntFieldsRelative)
      delete[] m_pnIntFieldsRelative;
  if (NULL != m_pnFloatFieldsRelative)
      delete[] m_pnFloatFieldsRelative;
  if (NULL != m_pnCstringFieldsRelative)
      delete[] m_pnCstringFieldsRelative;

  m_pcSelect            = NULL;
  m_nIntReflnFields     = 0;
  m_nFloatReflnFields   = 0;
  m_nCstringReflnFields = 0;
  m_nTotalFieldsPlus1   = 0;
}

int
Creflnlist::nInitValues(void)
{
  m_eThe_State         = eReflnlist_unknown_state;
  m_sReflnlist_key     = "unknown";

  

  m_fPercentRefToDisplay = 1.0f;
  m_bNeedsSaving         = 0;
  m_bOverwrite           = FALSE;
  m_eWriteBinary         = eReflnoutput_default;
  m_eWriteLog            = eLogoutput_default;
  m_nVerbose           = 1;
  m_nNumReflns         = 0;
  m_nNumReflnsAvail    = 0;
  m_nNumAlloc          = 0;
  m_nNumFree           = 0;
  m_nNumAllocFree      = 0;
  m_nAllocSize         = ALLOC_SIZE;
  m_poReflnMin         = NULL;
  m_poReflnMax         = NULL;
  m_ppoTheReflns       = NULL;
  m_ppoFreeReflns      = NULL;
  m_pnSortIndex        = NULL;
  m_poBinaryHeader     = NULL;
  m_nIntReflnFields    = 3;
  m_nFloatReflnFields  = 2;
  m_nCstringReflnFields= 0;
  m_nTotalFieldsPlus1  = m_nIntReflnFields + m_nFloatReflnFields
                                           + m_nCstringReflnFields + 1;
  m_psIntNames         = new Cstring [   m_nIntReflnFields ];
  m_psFloatNames       = new Cstring [ m_nFloatReflnFields ];
  if (0 < m_nCstringReflnFields)
    m_psCstringNames   = new Cstring [ m_nCstringReflnFields ];
  else
    m_psCstringNames   = NULL;

  m_pcSelect           = new char [m_nTotalFieldsPlus1];
  for (int i = 0; i < m_nTotalFieldsPlus1; i++) *(m_pcSelect+i) = (char) 0;

  *m_psIntNames         = ms_snH;
  *(m_psIntNames+1)     = ms_snK;
  *(m_psIntNames+2)     = ms_snL;
  *m_psFloatNames       = ms_sfIntensity;
  *(m_psFloatNames+1)   = ms_sfSigmaI;

  m_pnIntFieldsRelative = NULL;
  m_pnFloatFieldsRelative = NULL;
  m_pnCstringFieldsRelative = NULL;

  m_nNoWriteHeader      = 0;
  m_pnFortranFormat     = NULL;
  m_pnFortranFormat2    = NULL;

  vUpdateFieldIndices();
  nEraseExtra();

  return 0;
}

int
Creflnlist::nList(const int nNum)
{
  int nNumToList;
  Crefln *poRefln;
  cout << "Creflnlist::nList called.\n";
  cout << "     There are " << m_nNumReflns << " reflections in the list.\n";

  nNumToList = min(nNum, m_nNumReflns);
  for (int i = 0; i < nNumToList; i++)
    {
      poRefln = *(m_ppoTheReflns+i);
      cout << "Refln " << i << ": "
           << poRefln->nGetH() << " "
           << poRefln->nGetK() << " "
           << poRefln->nGetL() << endl;
    }
  return (0);
}

int
Creflnlist::nListFields(void)
{
  int i;
  cout << " Names of reflection fields:\n";
  for (i = 0; i < m_nIntReflnFields; i++)
    cout << "   " << *(m_psIntNames+i) << " (int)\n";
  for (i = 0; i < m_nFloatReflnFields; i++)
    cout << "   " << *(m_psFloatNames+i) << " (float)\n";
  for (i = 0; i < m_nCstringReflnFields; i++)
    cout << "   " << *(m_psCstringNames+i) << " (Cstring)" << endl;
  return (0);
}


int Creflnlist::nGetNumReflns(const Cstring& sFilename) {
    Creflnlist oList;
    
    oList.m_nNumReflnsAvail = -1;
    oList.nRead(sFilename,1,0);
    return oList.m_nNumReflnsAvail;
}

int
Creflnlist::nRead(const Cstring& sFilename, int nMaxToRead, int nFirstToRead, bool bIgnoreReadError)
{
    int                           nStat;
    int                           nx,ny;
    Cstring                       sTemp;
    int                           nReadIn;
    FILE*                         pFIn = NULL;               // Input file stream.
    char*                         pcHeadData = NULL;         // Header data read into this buffer.
    int                           nHeadData;          // Number of bytes in header data buffer.
    tagBinary_Ref_Header          oRefHeader;         // Binary reflection header (if used).
    char*                         pcExtra;            // Extra data (if used).
    
    if (0 != m_nNumReflns)
    {
        cout << "Creflnlist::nRead must read into an EMPTY Creflnlist!\n";
        return (1);
    }
    
    if (0 < m_nVerbose)
        cout << "Creflnlist::nRead with filename: " << sFilename << '\n';
    
    sTemp = sTransSymbol(sFilename);
    pFIn = fopen(sTemp.string(),"rb");
    
    if( NULL == pFIn )
    {
        cout << "ERROR: could not open " << sFilename << '\n';
		cout << strerror(errno) << std::endl; //mrp
        return 1;
    }
    
    // Read the header data.      
    // If the header data exceeds 10240 in size, then there may be a bug
    // on Mac OSX systems
    
    // RB 07/05/07 Actually, if the header data exceeds 10240 in size, there will be a problem 
    // on ANY operating system. The fix is to re-allocate the buffer, once 
    // the algorithm detects that 10240 bytes is not enough. That is done in 
    // tagBinary_Ref_Header::nReadHeader()
    
    nHeadData = max(200, min(10240, lFileGetSize(sTemp)));
    pcHeadData = new char[nHeadData];
    memset(pcHeadData, 0, nHeadData);
    nStat = fread(pcHeadData, nHeadData-1, 1, pFIn);
    fclose(pFIn);
    
    nReadIn = 0;
    
    if( nStat == 1 )
    {             
        nStat = nReadHeader(sFilename,
                            &pcHeadData,
                            nHeadData,
                            &pcExtra,
                            oRefHeader,
                            m_nBytesOrLinesReadInHeader);
        
        // If nStat==1 Could not read the text header.
        // If nStat==-1 then we failed to find the ID string. (Thus not a binary file)
        // If nStat==0 then we successfully read the header and did an integrity check.      
        
        if( nStat > 0 )
        {
            cout << "Could not read text header for file " << sFilename << "\n" << flush;
            delete[] pcHeadData;
            return  nStat;
        }
        
        if (1 < m_nVerbose)
            (void) nListFields(); // debuggin
        
        // Now create dummy reflection with all necessary fields.  If it had been
        // created previously, then 'this' list would not have had all the
        // necessary fields.
        
        Crefln oRefln(this);
        
        if (! nStat) {
            
            // Read in the extra data.
            if (pcExtra) {
                nParseExtra(pcExtra);
                delete[] pcExtra;
            }
            
            pFIn = fopen(sTransSymbol(sFilename).string(),"rb");
            nStat = 1;
            
            if (pFIn) {
                nStat = 0;
                
                for (nx=0;nx<m_nBytesOrLinesReadInHeader;nx++) {
                    if (1!=fread(&ny,1,1,pFIn))
                        nStat++;
                }
                
                m_eWriteBinary = eReflnoutput_binary; 
                
                
                if ((nStat==0) && (nMaxToRead != 0)) {
                    // Find how much data there remains to be read.
                    int nPosInit;
                    int nPosEnd;
                    int nSize;
                    int nMaxReflectionWidth;
                    int nMemorySize;
                    int nMemoryPointer;
                    int nMemoryReadIn;
                    int nMemoryProcessed;
                    char* pcBuf;
                    char* pcBufPoint;
                    
                    nPosInit = ftell(pFIn);
                    fseek(pFIn,0,SEEK_END);
                    nPosEnd = ftell(pFIn);
                    fseek(pFIn,nPosInit,SEEK_SET);
                    nSize = nPosEnd - nPosInit;
                    nMaxReflectionWidth = ((int) oRefHeader.cSizeInt)*oRefHeader.a3nFields[0]+((int) oRefHeader.cSizeFloat)*oRefHeader.a3nFields[1]+oRefHeader.a3nFields[2]*oRefHeader.nMaxString;
                    
                    if (nMaxToRead==-1) {
                        // We place a little padding at the end.  In the case that we read in bogus data,
                        // this will hopefully keep the system from core dumping.
                        nMemorySize = nSize + nMaxReflectionWidth;
                    } else {
                        // We will only be reading a portion of the data.
                        nMemorySize = nMaxReflectionWidth*nMaxToRead;
                    }
                    nMemorySize = min(1000000,nMemorySize);
                    pcBuf = new char[nMemorySize + 10];                 
                    
                    nReadIn = 0;                  
                    pcBufPoint = pcBuf;
                    nMemoryPointer = nMemorySize;
                    nMemoryReadIn = 0;
                    nMemoryProcessed = 0;
                    
                    while (((nMaxToRead==-1) || (nReadIn<nMaxToRead + nFirstToRead)) && (nMemoryProcessed<nSize)) {
                        
                        if (nMemoryPointer >= nMemorySize - nMaxReflectionWidth*3) {
                            if (nMemoryReadIn<nSize) {
                                // Move the memory.
                                memmove(pcBuf,pcBufPoint,nMemorySize - nMemoryPointer);
                                
                                // Read in up to the limit of the memory buffer.
                                nx =  fread(pcBuf + nMemorySize - nMemoryPointer,1,nMemoryPointer,pFIn);
                                if (nx == 0) {
                                    printf("ERROR:  Could not read in all data.  Had to terminate after %d reflections!\n",nReadIn);
                                    break;
                                }
                                nMemoryReadIn += nx;
                                
                                // Reset variables.
                                nMemoryPointer = 0;
                                pcBufPoint = pcBuf;
                            }
                        }
                        
                        nx = oRefHeader.nRead(pcBufPoint,oRefln);
                        nMemoryProcessed += nx;
                        nMemoryPointer += nx;
                        pcBufPoint += nx;
                        
                        if (nReadIn>=nFirstToRead) 
                            nStat = nInsert(&oRefln);
                        else
                            nStat = 0;

                        if (nStat != 0) {
                            printf("ERROR in Creflnlist::nRead inserting.\n");
                            break;
                        }
                        nReadIn++;
                    }
                    
                    if ((nMemoryProcessed != nSize) && (nReadIn!=nMaxToRead + nFirstToRead)) {
                        printf("ERROR:  Could not read entire reflection list!\n");
//+JWP 2007-09-07
			// We sometimes want to ignore this error
			if (bIgnoreReadError)
			  {
			    printf("INFO in Creflnlist::nRead read error ignored.\n");
			    nStat = 0;
			  }
			else
			  nStat = 1;
//-JWP 2007-09-07
                    } else
                        nStat = 0;

                    nReadIn = nGetNumReflns();
                    
                    if (m_nNumReflnsAvail == -1)
                        m_nNumReflnsAvail = oRefHeader.nNumRefs;

                    if (oRefHeader.nNumRefs==0) {
                        oRefHeader.nNumRefs = nReadIn;                        
                    } else 
                        nStat = -1;
                    if (0 < m_nVerbose) {
                        cout << "INFO in Creflnlist::nRead, EOF after "
                            << nReadIn << " reflections read in\n"
                            << "                                    ("
                            << m_nNumReflns
                            << " total now in list).\n";
                    }
                    
                    
                    delete[] pcBuf;
                }
          }
          
     } else if (nStat == -1) {
         m_eWriteBinary = eReflnoutput_text; 
         nStat = 0;
         int nCount=0;
         
         pFIn = fopen(sTransSymbol(sFilename).string(),"rt");
         nStat = 1;
         if (pFIn) {
             nStat = 0;
             
             for (nx=0;nx<m_nBytesOrLinesReadInHeader;nx++) {
                 if ((!fgets(pcHeadData,nHeadData,pFIn)) ||
                     (pcHeadData[strlen(pcHeadData)-1]!='\n'))
                     nStat++;
             }

             while (0 == nStat && ((nMaxToRead == -1) || (nCount < nMaxToRead + nFirstToRead)))
             {
                 if (nCount < nFirstToRead) {
                     char pcBuf[1000];
                     if ((NULL == fgets(pcBuf,1000-1,pFIn)) || (pcBuf[strlen(pcBuf)-1]!='\n'))
                         nStat = 1;
                     else
                         nStat = 0;

                 } else
                     nStat = oRefln.nRead(pFIn);

                 if (0 == nStat)
                 {
                     if (nCount >= nFirstToRead) {
                         nReadIn++;
                         nStat = nInsert(&oRefln);
                         if (nStat != 0)
                         {
                             cout << "ERROR in Creflnlist::nRead inserting." << endl;
                         }
                     }
                 }
                 else if (feof(pFIn))
                 {
                     if (0 < m_nVerbose)
                     {
                         cout << "INFO in Creflnlist::nRead, EOF after "
                             << nReadIn << " reflections read in\n"
                             << "                                    ("
                             << m_nNumReflns
                             << " total now in list).\n";
                     }
                     nStat = -1;
                 }
                 else
                 {
                     cout << "ERROR in Creflnlist::nRead, error after "
                         << nReadIn << " reflections read in (" << m_nNumReflns
                         << " total now in list).\n";
                     nStat = 3;
                 }
                 nCount++;
             }
             if ((m_nNumReflnsAvail == -1) && ((nStat==-1) || (nStat==0))) {
                 nStat = 0;
                 nCount = 0;
                 while ((!feof(pFIn)) && (!nStat)) {
                     char pcBuf[1000];
                     if ((NULL == fgets(pcBuf,1000-1,pFIn)) || (pcBuf[strlen(pcBuf)-1]!='\n'))
                         nStat = 1;
                     else {
                         nStat = 0;
                         nCount++;
                     }                     
                 }
                 m_nNumReflnsAvail = nCount;
             }


         }
     }
     fclose(pFIn);
     
     
     if (-1 == nStat) 
         nStat = 0;
    }
    else
    {
        cout << "ERROR in Creflnlist::nRead opening: " << sFilename << '\n';
        nStat = 1;
    }
    vUpdateFieldIndices();
    
    if(m_nFI_fObsPx0 > 0)
    {
        m_nObservedSpotsExist = 1;
    }
    else
    {
        m_nObservedSpotsExist = 0;
    }
    if(m_nFI_fCalcPx0 > 0)
    {
        m_nCalculatedSpotsExist = 1;
    }
    else
    {
        m_nCalculatedSpotsExist = 0;
    }
    
    if (NULL != pcHeadData) delete[] pcHeadData;
    return (nStat);
}


int
Creflnlist::nWrite(const Cstring& sFilename, int *pnIndex,
                   const char *pcSelectIn,
                   FILE* pFOutUser,const int nFirst, const int nLast)
{
  Crefln *poRefln;
  int i;
  int *pnTemp;
  int nStat;
  int nFirstLocal;
  int nLastLocal;
  int nNumWritten = 0;
  FILE* pFOut;
  static Cstring sTemp;


  // If the user specified an output stream, then we should
  // ignore the input filename.
  
  if (pFOutUser != NULL ) {
      pFOut = pFOutUser;
      
  } else {
      
    // If the member variable for overwriting is not set, 
    // then check if the overwrite environment variable is set, then make sure
    // any existing file gets a version number appended, so it is
    // not overwritten
      
    if (m_bOverwrite)
      nStat = 0;
    else
      nStat = nFileAppendVersion(sFilename, TRUE);
      
    if (0 != nStat)
      {
          // Reflection lists are so important, this is a FATAL error
          
          cout << "FATAL ERROR renaming reflnlist file: " << sFilename << "!\n";
#ifdef SSI_PC
          return nStat;
#else
          exit (nStat);
#endif
      }
      // Fortran format will differ depending upon the machine on which it is generated.
      // Really, we should have a seperate flag to control this "if", but it works for now.
      if (m_pnFortranFormat)
          pFOut = fopen(sTransSymbol(sFilename).string(),"wt");
      else
          pFOut = fopen(sTransSymbol(sFilename).string(),"wb");
  }


  bool bNewIndex = FALSE;

  if (NULL == pnIndex)
    {
      // No index array for the order of the reflections was provided, so
      // use internal one if not NULL, otherwise create one.

      if (NULL == m_pnSortIndex)
        {
          pnTemp = new int [m_nNumReflns];
          bNewIndex = TRUE;
          for (i = 0; i < m_nNumReflns; i++) pnTemp[i] = i;
        }
      else
        {
          pnTemp = m_pnSortIndex;
        }
    }
  else
    {
      pnTemp = pnIndex;
    }

  nFirstLocal = 0;
  nLastLocal  = m_nNumReflns-1;
  if (0 <= nFirst)
    {
      nFirstLocal = nFirst;
      nFirstLocal = min(nFirstLocal, m_nNumReflns-1);
    }
  if (0 <= nLast)
    {
      nLastLocal = nLast;
      nLastLocal = min(nLastLocal, m_nNumReflns-1);
    }

  // Delete the old binary header if it exists.
  if (NULL != m_poBinaryHeader) 
      delete m_poBinaryHeader;
  m_poBinaryHeader=NULL;

  nStat = 0;
  if (pFOut)
    {
      bool bWriteBinary;

      if (m_eWriteBinary==eReflnoutput_binary)
          bWriteBinary=TRUE;
      else if (m_eWriteBinary==eReflnoutput_text)
          bWriteBinary=FALSE;
      else if (sGetEnv("DTREK_REFLN_BINARY").length()!=0) {
          bWriteBinary=TRUE;          
      } else 
          bWriteBinary=FALSE;

      // Write the TEXT header.  This goes in binary and text style formats.
      nStat = nWriteHeader(pFOut, pcSelectIn,bWriteBinary);

      if (nStat) {
            cout << "ERROR writing header in Creflnlist::nWrite\n" << flush;
      }

      if ((nStat == 0) && (bWriteBinary)) {
          // Create a new binary header.  
          m_poBinaryHeader = new  tagBinary_Ref_Header;
          // Write the binary header fields.
          m_poBinaryHeader->nLoadHeader(*this,pcSelectIn);
          // Fix up the number of reflections if we are writing only a portion of the list.
          if ((nFirst!=-1) && (nLast!=-1))
            m_poBinaryHeader->nNumRefs = nLast - nFirst + 1;
          // We need the size of the extra data.
          nGenerateLog(sTemp,sFilename,*m_poBinaryHeader);
          m_poBinaryHeader->nExtraSize = sTemp.length()+1+m_sEmbeddedHeader.length()+1;
          m_poBinaryHeader->nWriteHeader(pFOut);
          // Write the extra data.
          fwrite(sTemp.string(),sTemp.length()+1,1,pFOut);
          fwrite(m_sEmbeddedHeader.string(),m_sEmbeddedHeader.length()+1,1,pFOut);
      }
      
      // Lines 3-?: Write out selected reflections and selected fields
      
      if (nFirstLocal != -1) {
          for (i = nFirstLocal; (i <= nLastLocal) && (0 == nStat); i++)
          {
              poRefln = m_ppoTheReflns[pnTemp[i]];
              poRefln->vWrite(pFOut,pcSelectIn);
              
              if (ferror(pFOut))
              {
                  cout << "ERROR in Creflnlist::nWrite writing: "
                      << sFilename << cEOL; 
                  nStat = 1;
              }
              else
              {
                  nNumWritten++;
              }
          } // End of for... loop
      }
      
      // If the stream was not provided on the command line, then we need to shut it.
      // Otherwise, the user has more reflection fields to write.
      if (pFOutUser == NULL) {
#ifdef WIN32
		    int iFileDes = fileno(pFOut);  //Added KWY 10/7/09
			_commit(iFileDes);			   //Added KWY 10/7/09
#endif
            fclose(pFOut);
#ifdef WIN32
			Cstring csDirFile = sFilename.sGetFilePath();  //Added KWY 10/20/09
			csDirFile += "\\.";						       //Added KWY 10/20/09
			FILE* fpDirFile = fopen(csDirFile, "r+");	   //Added KWY 10/20/09
			if( fpDirFile ) {							   //Added KWY 10/20/09
				int iDirDes = fileno(fpDirFile);		   //Added KWY 10/20/09
				_commit(iDirDes);						   //Added KWY 10/20/09
				fclose(fpDirFile);						   //Added KWY 10/20/09
			}
#endif
      }
    }
  else
    {
      cout << "ERROR in Creflnlist::nWrite opening: " << sFilename << cEOL;
      nStat = 2;
    }
  if (bNewIndex)
    {
      delete [] pnTemp;
      pnTemp = NULL;
    }
  m_bNeedsSaving = 0;

  if ( (pFOutUser != NULL ) && (0 == nNumWritten) )
  {
      printf("... reflnlist %s prepared for writing.\n",sFilename.string());
  }
  else
      printf("Number of reflections written in '%s': %d\n",sFilename.string(),nNumWritten);

  cout << flush;  

  return (nStat);
}

int Creflnlist::nInsert(Crefln* poRefln, int nIndex)
{
  int nStat;

  if (nIndex == -1)
      nIndex = m_nNumReflns;
  else if ((nIndex<0) || (nIndex>m_nNumReflns))
      return 1;

  nStat = nExpand((const int) 1);

  if (0 == nStat)
    {
      // Move all of the existing reflections so that
      // the new one can be inserted at the specified index.
      if (nIndex < m_nNumReflns)
        memmove((void*) (&m_ppoTheReflns[nIndex+1]),(void*) (&m_ppoTheReflns[nIndex]),(m_nNumReflns - nIndex)*sizeof(Crefln*));
      
      // Integrate uses this copy a lot.  I replaced it with memmove to speed up processing on possibly large lists. (tjn)
      // for(i=m_nNumReflns; i>nIndex; --i)
      //   {
      //     m_ppoTheReflns[i] = m_ppoTheReflns[i-1];
      //  }

      // Add to list then.
      if (0 >= m_nNumFree)
        {
          // None available in free list, allocate new one
          m_ppoTheReflns[nIndex] = new Crefln (*poRefln);        // Danger here!
          m_ppoTheReflns[nIndex]->m_poReflnlist = this;
          // Because they most point back to same Reflnlist
          // for this to be valid?
        }
      else
        {
          // Refln available in free list, use it

          m_ppoTheReflns[nIndex] = m_ppoFreeReflns[m_nNumFree-1];
          *(m_ppoTheReflns[nIndex]) = *poRefln;
          m_ppoTheReflns[nIndex]->m_poReflnlist = this;
          m_ppoFreeReflns[m_nNumFree-1] = NULL;
          m_nNumFree--;
        }
      m_nNumReflns++;
    }
  return (nStat);
}


int
Creflnlist::nInsert(const int nH, const int nK, const int nL)
{
  int nStat;

  nStat = nExpand((const int) 1);

  if (0 == nStat)
    {
      // Add to list then.
      //   Perhaps get a poRefln object from a "free list"

      if (0 >= m_nNumFree)
        {
          // None available in free list, allocate new one

          m_ppoTheReflns[m_nNumReflns] = new Crefln(this, nH, nK, nL);
        }
      else
        {
          // Refln available in free list, use it
          // (Maybe set values to -999?

          m_ppoTheReflns[m_nNumReflns] = m_ppoFreeReflns[m_nNumFree-1];
          m_ppoTheReflns[m_nNumReflns]->vSetH(nH);
          m_ppoTheReflns[m_nNumReflns]->vSetK(nK);
          m_ppoTheReflns[m_nNumReflns]->vSetL(nL);
          m_ppoFreeReflns[m_nNumFree-1] = NULL;
          m_nNumFree--;
        }

      m_nNumReflns++;
    }
  return (nStat);
}

int  
Creflnlist::nInsertSorted(Crefln* poRefln,const eReflnFieldType eType,const int nField) {
    int nStat;
    int nInsertSpot;

    // We do NOT want to have a sorting index.  We assume that someone using this routine will 
    // want the reflections to occur in order WITHOUT reference to a sorting list.
    if (m_pnSortIndex)
        delete[] m_pnSortIndex;
    m_pnSortIndex = NULL;

    switch (eType) {
    case eReflnField_int_type:
        nInsertSpot = nFindFirst(nField,poRefln->nGetField(nField),NULL,TRUE);
        break;
    case eReflnField_float_type:
        nInsertSpot = nFindFirst(nField,(const double) poRefln->fGetField(nField),0.0,NULL,TRUE);
        break;
    case eReflnField_Cstring_type:
        nInsertSpot = nFindFirst(nField,poRefln->sGetField(nField),NULL,TRUE);
        break;
    }
    if (nInsertSpot == -1)
        nInsertSpot = m_nNumReflns;

    nStat = nInsert(poRefln,nInsertSpot);
    return (nStat);
}

int
Creflnlist::nExpand(const int nAdditional)
{
  Crefln **ppoNewList = NULL;
  int    i;
	  
  if (m_nNumAlloc <= (m_nNumReflns + nAdditional) )
    {

      // Need to allocate more memory, so create new list ...

      i = max(m_nNumAlloc + m_nAllocSize, m_nNumReflns + nAdditional);
      //+jwp 7-May-2001
      // Also expand by at least 10%
      i = max(i, m_nNumReflns + m_nNumReflns / 10);
      //-jwp 7-May-2001

      /******
      //+2012-May-22 TEST CODE

      ppoNewList = new Crefln* [i];
      cout << "BEFORE first delete, i is " << i << endl;
      delete [] ppoNewList;
      ppoNewList = NULL;
      cout << "AFTER first delete\n" << endl << flush;

      //-2012-May-22 TEST CODE
      ******/


      ppoNewList = new Crefln* [i];
      m_nNumAlloc = i;

      // Transfer Crefln pointers in old list to new list

      for (i = 0; i < m_nNumReflns; i++)
        {
          //*(ppoNewList + i) = *(m_ppoTheReflns + i);
          ppoNewList[i] = m_ppoTheReflns[i];
        }

      // Make sure the newly allocated ones point to NULL
      // since they do not exist

      for (i = m_nNumReflns; i < m_nNumAlloc; i++)
        ppoNewList[i] = NULL;

      // Clear old memory holding reflection list
          if (NULL != m_ppoTheReflns)
	    {
	      delete [] m_ppoTheReflns;
	      //+2012-05-22 jwp
	      m_ppoTheReflns = NULL;
	      //-2012-05-22 jwp
	    }

      //  ... and finally point to new memory location.

      m_ppoTheReflns = ppoNewList;
    }
  return (0);
}

int
Creflnlist::nSelectField(const eReflnFieldType eType,
                         const int nSelectField, const char cValue,
                         char *pcSelectIn)
{
  // Set the select character in pcSelectIn to the value specified in cValue
  //   for the specified (by eType and nSelectField) field.
  //   If pcSelectIn is NULL, use the internal member m_pcSelect variable.

  char *pcSelect;

  if ( (eReflnField_unknown_type == eType) || (0 > nSelectField) )
    {
      return (-1);  // Too hard to figure out, so return with error.
    }

  if (NULL == pcSelectIn)
    {
      // Use the internal selection character array since none was passed

      pcSelect = m_pcSelect;
    }
  else
    {
      pcSelect = pcSelectIn;
    }

  if (   (eReflnField_int_type == eType)
      && (nSelectField < m_nIntReflnFields) )
    {
      pcSelect[nSelectField] = cValue;
    }
  else if (   (eReflnField_float_type == eType)
           && (nSelectField < m_nFloatReflnFields) )
    {
      pcSelect[nSelectField + m_nIntReflnFields] = cValue;
    }
  else if (   (eReflnField_Cstring_type == eType)
           && (nSelectField < m_nCstringReflnFields) )
    {
      pcSelect[nSelectField + m_nIntReflnFields + m_nFloatReflnFields] = cValue;
    }
  else
    {
      return (-2);  // Some kind of error
    }
  return (0);
}

int
Creflnlist::nSelect(const int nFI, const Cstring& sOperation,
                    const float fTestValue, const bool bIncExc)
{
  int     i;
  Crefln *poRefln;
  char    cTemp;
  float   fTemp;
  int     nFI1, nFI2;
  int     nASMD = 0;           // Normal, Add, Subtract, Multiply, Divide mode

  // If nFI >= 1000000, that means two fields are encoded into nFI

  nASMD = nFI / 1000000;
  nFI1 = nFI % 1000;

  if (0 < nASMD)
    {
      nFI1 = nFI % 1000;
      nFI2 = (nFI - (nASMD * 1000000) - nFI1) / 1000;
      if ( (0 > nFI2) || (m_nFloatReflnFields <= nFI2) )
        return (-1); // Out-of-bounds
    }

  if ( (0 > nFI1) || (m_nFloatReflnFields <= nFI1) )
    return (-1); // Out-of-bounds

  y_prbLOGIC_F_F prbLogic = &bEqual;

  if (sOperation == ms_sOpEqual)
    {
      prbLogic = &bEqual;
    }
  else if (sOperation == ms_sOpNotEqual)
    {
      prbLogic = &bNotEqual;
    }
  else if (sOperation == ms_sOpLess)
    {
      prbLogic = &bLess;
    }
  else if (sOperation == ms_sOpGreater)
    {
      prbLogic = &bGreater;
    }
  else if (sOperation == ms_sOpGreaterEqual)
    {
      prbLogic = &bGreaterEqual;
    }
  else if (sOperation == ms_sOpLessEqual)
    {
      prbLogic = &bLessEqual;
    }
  else
    {
      // Illegal operation
      return (-1);
    }

  if (bIncExc)
    {
      // Include those where comparison is TRUE
      cTemp = char(0);
    }
  else
    {
      // Exclude those where comparison is TRUE
      cTemp = m_pcSelect[m_nTotalFieldsPlus1-1] + 1;
    }

  int nNumTrue = 0;
  static int s0DivideCount = 0;

  for (i = 0; i < m_nNumReflns; i++)
    {
      poRefln = *(m_ppoTheReflns+i);

      if (0 == nASMD)
        fTemp = poRefln->fGetField(nFI1);
      else if (1 == nASMD)
        fTemp = poRefln->fGetField(nFI1) + poRefln->fGetField(nFI2);
      else if (2 == nASMD)
        fTemp = poRefln->fGetField(nFI1) - poRefln->fGetField(nFI2);
      else if (3 == nASMD)
        fTemp = poRefln->fGetField(nFI1) * poRefln->fGetField(nFI2);
      else if (4 == nASMD)
        {
          // TODO: In order to prevent divide by zero, change divide to multiply

          fTemp = poRefln->fGetField(nFI2);
          if (0.0 != fTemp)
            {
              fTemp = poRefln->fGetField(nFI1) / fTemp;
            }
          else
            {
              s0DivideCount++;
              if (100 == s0DivideCount)
                {
                  cout << "WARNING in Creflnlist::nSelect: divisor is 0"
                       << " for at least 100 reflns!\n"
                       << flush;
                }

              // Often this is due to a (fIntensity/fSigmaI<fX) comparison,
              // so change to (-0.0001 < fX). 
              // Example: a -fIntensity/fSigmaI<5 will still exclude the
              //          a refln with fSigmaI of 0.

              fTemp      = -0.0001f;
            }
        }

      // Note that we test this reflection even if it was
      // was excluded by something else.

      if ((*prbLogic)(fTemp, fTestValue ))
        {
          // Include/Exclude this reflection

          poRefln->m_nSelect = cTemp;
          nNumTrue++;
        }
    }
  return (nNumTrue);
}

int
Creflnlist::nSelect(const int nFI, const Cstring& sOperation,
                    const int nTestValue, const bool bIncExc)
{
  // This actually excludes reflections which do not match the
  // criteria.  Excluding reflections is cumulative, that is, once
  // flagged excluded there is no code to unexclude them (for now).

  int     i;
  Crefln *poRefln;
  char    cTemp;
  int     nTemp;
  int     nFI1, nFI2;
  int     nASMD = 0;           // Normal, Add, Subtract, Multiply, Divide mode

  // If nFI >= 1000000, that means two fields are encoded into nFI

  nASMD = nFI / 1000000;
  nFI1 = nFI % 1000;

  if (0 < nASMD)
    {
      nFI1 = nFI % 1000;
      nFI2 = (nFI - (nASMD * 1000000) - nFI1) / 1000;
      if ( (0 > nFI2) || (m_nIntReflnFields <= nFI2) )
        return (-1); // Out-of-bounds
    }

  if ( (0 > nFI1) || (m_nIntReflnFields <= nFI1) )
    return (-1); // Out-of-bounds

  y_prbLOGIC_N_N  prbLogic = &bEqual;

  if (sOperation == ms_sOpEqual)
    {
      prbLogic = &bEqual;
    }
  else if (sOperation == ms_sOpNotEqual)
    {
      prbLogic = &bNotEqual;
    }
  else if (sOperation == ms_sOpLess)
    {
      prbLogic = &bLess;
    }
  else if (sOperation == ms_sOpGreater)
    {
      prbLogic = &bGreater;
    }
  else if (sOperation == ms_sOpGreaterEqual)
    {
      prbLogic = &bGreaterEqual;
    }
  else if (sOperation == ms_sOpLessEqual)
    {
      prbLogic = &bLessEqual;
    }
  else if (sOperation == ms_sOpAnd)
    {
      prbLogic = &bAnd;
    }
  else if (sOperation == ms_sOpOr)
    {
      prbLogic = &bOr;
    }
  else
    {
      // Illegal operation
      return (-1);
    }
  if (bIncExc)
    {
      // Include those where comparison is TRUE
      cTemp = char(0);
    }
  else
    {
      // Exclude those where comparison is TRUE
      cTemp = m_pcSelect[m_nTotalFieldsPlus1-1] + 1;
    }

  int nNumTrue = 0;
  for (i = 0; i < m_nNumReflns; i++)
    {
      poRefln = *(m_ppoTheReflns+i);

      if (0 == nASMD)
        nTemp = poRefln->nGetField(nFI1);
      else if (1 == nASMD)
        nTemp = poRefln->nGetField(nFI1) + poRefln->nGetField(nFI2);
      else if (2 == nASMD)
        nTemp = poRefln->nGetField(nFI1) - poRefln->nGetField(nFI2);
      else if (3 == nASMD)
        nTemp = poRefln->nGetField(nFI1) * poRefln->nGetField(nFI2);
      else if (4 == nASMD)
        nTemp = poRefln->nGetField(nFI1) / poRefln->nGetField(nFI2);

      // Note that we test this reflection even if it was
      // was excluded by something else.

      if ((*prbLogic)(nTemp, nTestValue ))
        {
          // Include/Exclude this reflection
          poRefln->m_nSelect = cTemp;
          nNumTrue++;
        }
    }
  return (nNumTrue);
}

int
Creflnlist::nSelect(const int nFI, const Cstring& sOperation,
                    const Cstring& sTestValue, const bool bIncExc)
{
  int i;
  Crefln *poRefln;
  char    cTemp;

  if ( (0 > nFI) || (m_nCstringReflnFields <= nFI) )
    return (-1); // Out-of-bounds

  y_prbLOGIC_S_S  prbLogic = &bEqual;

  if (sOperation == ms_sOpEqual)
    {
      prbLogic = &bEqualRegular;
    }
  else if (sOperation == ms_sOpNotEqual)
    {
      prbLogic = &bNotEqualRegular;
    }
  else if (sOperation == ms_sOpLess)
    {
      prbLogic = &bLess;
    }
  else if (sOperation == ms_sOpGreater)
    {
      prbLogic = &bGreater;
    }
  else if (sOperation == ms_sOpGreaterEqual)
    {
      prbLogic = &bGreaterEqual;
    }
  else if (sOperation == ms_sOpLessEqual)
    {
      prbLogic = &bLessEqual;
    }
  else
    {
      // Illegal operation
      return (-1);
    }
  if (bIncExc)
    {
      // Include those where comparison is TRUE
      cTemp = char(0);
    }
  else
    {
      // Exclude those where comparison is TRUE
      cTemp = m_pcSelect[m_nTotalFieldsPlus1-1] + 1;
    }

  int nNumTrue = 0;
  int nLen;
  nLen = sTestValue.length() - 1;
  if ('*' != sTestValue.GetAt(nLen))
    {
      // The last character in sTestValue is not an asterisk (wildcard),
      // so look for exact match

      for (i = 0; i < m_nNumReflns; i++)
        {
          poRefln = *(m_ppoTheReflns+i);
        
          // Note that we test this reflection even if it was
          // was excluded by something else.
        
          if ((*prbLogic)(poRefln->sGetField(nFI), sTestValue))
            {
              // Include/Exclude this reflection
              poRefln->m_nSelect = cTemp;
              nNumTrue++;
            }
        }
    }
  else
    {
      // The last character in sTestValue is an asterisk, so look for
      // matches up to that character.

      Cstring sTemp1, sTemp2;
      sTemp2 = Cstring(sTestValue, 0, nLen);

      for (i = 0; i < m_nNumReflns; i++)
        {
          poRefln = *(m_ppoTheReflns+i);
        
          // Note that we test this reflection even if it was
          // was excluded by something else.
        
          sTemp1 = Cstring(poRefln->sGetField(nFI), 0, nLen);
          if ((*prbLogic)(sTemp1, sTemp2))
            {
              // Include/Exclude this reflection
              poRefln->m_nSelect = cTemp;
              nNumTrue++;
            }
        }
    }
  return (nNumTrue);
}

int
Creflnlist::nSelect(const Cstring& sCompareString)
{
  // Include/Exclude reflections based on sCompareString of the
  // form  [+|-]fieldName1[MfieldName2]OPERATORvalue (no whitespace!)
  // where + means include
  //       - means exclude
  //       fieldName1 is a fieldname in the reflection list
  //       M          is a character +, -, *, / for a math operation
  //       fieldName2 is a fieldname in the reflection list of same type as fn1
  //       OPERATOR  is one of the operator strings of class Creflnlist
  //       value     is a valid value
  //
  // or the special sCompareString:  -rejectbatchid=template

  int      i;
  Cstring sTemp, sBefore, sBeforeB, sOp, sAfter;
  bool bIncExc;

  sTemp      = sCompareString;   // Get rest of string

  if ("-1to2" == sTemp)
    {
      return (nConvertFromAnom());
    }
  else if (sTemp.contains("-rejectbatchid="))
    {
      sTemp = sTemp.after("-rejectbatchid=");
      return (nDeleteBatchID(sTemp));
    }

  if ('+' == sTemp.GetAt(0))
    {
      bIncExc = TRUE;
    }
  else if ('-' == sTemp.GetAt(0))
    {
      bIncExc = FALSE;
    }
  else
    {
      cout << "Creflnlist::nSelect(" << sCompareString << ") is illegal!\n";
      return (-1);    // Problem parsing
    }

  sTemp = sTemp.after(0);  // Remove first character

  if (0 < sTemp.index(ms_sOpEqual))
    {
      // == found
      sOp = ms_sOpEqual;
    }
  else if (0 < sTemp.index(ms_sOpEqual_Fortran))
  {
      sOp = ms_sOpEqual_Fortran;
  }
  else if (0 < sTemp.index(ms_sOpNotEqual))
    {
      // != found
      sOp = ms_sOpNotEqual;
    }
  else if (0 < sTemp.index(ms_sOpNotEqual_Fortran))
    {
      // != found
      sOp = ms_sOpNotEqual_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpGreaterEqual))
    {
      // >= found
      sOp = ms_sOpGreaterEqual;
    }
  else if (0 < sTemp.index(ms_sOpGreaterEqual_Fortran))
    {
      // >= found
      sOp = ms_sOpGreaterEqual_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpLessEqual))
    {
      // <= found
      sOp = ms_sOpLessEqual;
    }
  else if (0 < sTemp.index(ms_sOpLessEqual_Fortran))
    {
      // <= found
      sOp = ms_sOpLessEqual_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpAnd))
    {
      // & found
      sOp = ms_sOpAnd;
    }
  else if (0 < sTemp.index(ms_sOpAnd_Fortran))
    {
      // & found
      sOp = ms_sOpAnd_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpOr))
    {
      // | found
      sOp = ms_sOpOr;
    }
  else if (0 < sTemp.index(ms_sOpOr_Fortran))
    {
      // | found
      sOp = ms_sOpOr_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpAddEqual))
    {
      // += found
      sOp = ms_sOpAddEqual;
    }
  else if (0 < sTemp.index(ms_sOpSubEqual))
    {
      // -= found
      sOp = ms_sOpSubEqual;
    }
  else if (0 < sTemp.index(ms_sOpMulEqual))
    {
      // *= found
      sOp = ms_sOpMulEqual;
    }
  else if (0 < sTemp.index(ms_sOpDivEqual))
    {
      // /= found
      sOp = ms_sOpDivEqual;
    }
  else if (0 < sTemp.index(ms_sOpSetEqual))
    {
      // = found
      sOp = ms_sOpSetEqual;
    }
  else if (0 < sTemp.index(ms_sOpLess))
    {
      // < found
      sOp = ms_sOpLess;
    }
  else if (0 < sTemp.index(ms_sOpLess_Fortran))
    {
      // < found
      sOp = ms_sOpLess_Fortran;
    }
  else if (0 < sTemp.index(ms_sOpGreater))
    {
      // > found
      sOp = ms_sOpGreater;
    }
  else if (0 < sTemp.index(ms_sOpGreater_Fortran))
    {
      // > found
      sOp = ms_sOpGreater_Fortran;
    }
  else if ("nonaxial" == sTemp)
    {
      // Special +-nonaxial selection

      int      nNumNotAxial = 0;
      int      nNumZeroes;
      Crefln *poRefln;
      char     cTemp;

      if (bIncExc)
        {
          // Include those where comparison is TRUE
          cTemp = char(0);
          cout << "WARNING +axes not supported!\n" << flush;
          return (-1);
        }
      else
        {
          // Exclude those where comparison is TRUE
          cTemp = m_pcSelect[m_nTotalFieldsPlus1-1] + 1;
        }
      for (i = 0; i < m_nNumReflns; i++)
        {
          poRefln = poGetRefln(i);
          nNumZeroes = 0;
          if (0 == poRefln->nGetH()) nNumZeroes++;
          if (0 == poRefln->nGetK()) nNumZeroes++;
          if (0 == poRefln->nGetL()) nNumZeroes++;
          if (2 > nNumZeroes)
            {
              // It is not an axial refln, exclude it

              poRefln->m_nSelect = cTemp;
              nNumNotAxial++;
            }
        }
      return (nNumNotAxial);
    }
  else
    {
      cout << "Comparison operator not found!\n";
      cout << "Creflnlist::nSelect(" << sCompareString << ") is illegal!\n";
      return (-1);
    }

  sBefore = sTemp.before(sOp);
  sAfter  = sTemp.after(sOp);

  //  Look for +, -, *, / in the fieldname which indicates 2 fieldnames given
  int  nASMD = 0;
  for (i = 0; (i < 4) && (0 >= nASMD); i++)
    {
      nASMD = sBefore.index(ms_sOpASMD.GetAt(i));
      if (0 < nASMD)
        {
          // Found a match, so separate out before and after
          sBeforeB = sBefore.after((int)nASMD);
          sBefore  = sBefore.before((int)nASMD);
          nASMD    = i + 1;
        }
    }

  // Ok, have before, operator, and after

  int      nTemp, nGot;
  int      nFI, nFI2;
  float    fTemp;
  eReflnFieldType eType = eReflnField_unknown_type;

  // What type field is it?  Use the first character of the
  // field name first.  If that does not work,
  // then use the value itself.

  nGot = 0;
  if ('n' == sBefore.GetAt(0))
    {
      eType = eReflnField_int_type;
//      nGot = sscanf((const char *)sAfter, "%d", &nTemp);
      nGot = sscanf(sAfter.string(), "%d", &nTemp);
    }
  else if ('f' == sBefore.GetAt(0))
    {
      eType = eReflnField_float_type;
//      nGot = sscanf((const char *)sAfter, "%f", &fTemp);
      nGot = sscanf(sAfter.string(),"%f", &fTemp);
    }
  else if ('s' == sBefore.GetAt(0))
    {
      eType = eReflnField_Cstring_type;
      if (0 < sAfter.length())
        {
          nGot = 1;
        }
    }
  else
    {
      // Try to figure type out from the value


//      nGot = sscanf((const char *)sAfter, "%d", &nTemp);
      nGot = sscanf(sAfter.string(), "%d", &nTemp);
      if ( (1 == nGot) && (0 > sAfter.index('.'))     // No decimal point
                       && (0 > sAfter.index('E'))     // No E as in:  1E10
                       && (0 > sAfter.index('e')) )   // No e as in:  1e10
        {
          eType = eReflnField_int_type;
        }
      else
        {
//          nGot = sscanf((const char *)sAfter, "%f", &fTemp);
          nGot = sscanf(sAfter.string(), "%f", &fTemp);
          if (1 == nGot)
            {
              eType = eReflnField_float_type;
            }
          else
            {
              eType = eReflnField_Cstring_type;
              nGot = 1;
            }
        }
    }
  if (1 != nGot)
    {
      cout << "Creflnlist::nSelect(" << sCompareString << ") is illegal!\n";
      return (-2);   // Another problem parsing!
    }

  nFI = nGetFieldIndex(sBefore, eType);
  if (0 < nASMD)
    {
      nFI2 = nGetFieldIndex(sBeforeB, eType);
      if (0 > nFI2)
        {
          cout << "Fieldname " << sBeforeB << " not found!\n";
          return (-3);    // Could not find second fieldname
        }
      else
        {
          nFI = (1000000 * nASMD) + (1000 * nFI2) + nFI;
        }
    }

  // Convert fortran operators, back to standard operators.
  if (sOp == ms_sOpEqual_Fortran)
      sOp = ms_sOpEqual;
  else if (sOp == ms_sOpLessEqual_Fortran)
      sOp = ms_sOpLessEqual;
  else if (sOp == ms_sOpGreaterEqual_Fortran)
      sOp = ms_sOpGreaterEqual;
  else if (sOp == ms_sOpNotEqual_Fortran)
      sOp = ms_sOpNotEqual;
  else if (sOp == ms_sOpLess_Fortran)
      sOp = ms_sOpLess;
  else if (sOp == ms_sOpGreater_Fortran)
      sOp = ms_sOpGreater;
  else if (sOp == ms_sOpAnd_Fortran)
      sOp = ms_sOpAnd;
  else if (sOp == ms_sOpOr_Fortran)
      sOp = ms_sOpOr;

  if ( (0 > nFI) && (ms_sOpSetEqual != sOp) )
    {
      cout << "Creflnlist::nSelect(" << sCompareString << ") is illegal!\n";
      return (-3); // Another problem parsing!
    }

  if (ms_sOpSetEqual == sOp)
    {
      nFI = nExpandGetField(sBefore, eType);  // Insert field if missing
      if (eReflnField_int_type == eType)
        {
          // Integer field
          for (i = 0; i < m_nNumReflns; i++)
          {
            poGetRefln(i)->vSetField(nFI, nTemp);
          }
        }
      else if (eReflnField_float_type == eType)
        {
          for (i = 0; i < m_nNumReflns; i++)
            {
              poGetRefln(i)->vSetField(nFI, fTemp);
            }
        }
      else if (eReflnField_Cstring_type == eType)
        {
          for (i = 0; i < m_nNumReflns; i++)
            {
              poGetRefln(i)->vSetField(nFI, sAfter);
            }
        }
      return (m_nNumReflns);
    }
  else if (   (ms_sOpAddEqual == sOp)
           || (ms_sOpSubEqual == sOp)
           || (ms_sOpMulEqual == sOp)
           || (ms_sOpDivEqual == sOp) )
    {
      int nOp = -1;                             // Starting to look kludgy!
      if      (ms_sOpAddEqual == sOp) nOp = 0;
      else if (ms_sOpSubEqual == sOp) nOp = 1;
      else if (ms_sOpMulEqual == sOp) nOp = 2;
      else if (ms_sOpDivEqual == sOp) nOp = 3;

      if (eReflnField_int_type == eType)
        {
          // Integer field

          int nValue;
          for (i = 0; i < m_nNumReflns; i++)
            {
              nValue = poGetRefln(i)->nGetField(nFI);
              if (0 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, nValue+nTemp);
                }
              else if (1 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, nValue-nTemp);
                }
              else if (2 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, nValue*nTemp);
                }
              else if (3 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, nValue/nTemp);
                }
            }
        }
      else if (eReflnField_float_type == eType)
        {
          float fValue;
          for (i = 0; i < m_nNumReflns; i++)
            {
              fValue = poGetRefln(i)->fGetField(nFI);
              if (0 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, fValue+fTemp);
                }
              else if (1 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, fValue-fTemp);
                }
              else if (2 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, fValue*fTemp);
                }
              else if (3 == nOp)
                {
                  poGetRefln(i)->vSetField(nFI, fValue/fTemp);
                }
            }
        }
      else if (eReflnField_Cstring_type == eType)
        {
          if (0 == nOp)
            {
              for (i = 0; i < m_nNumReflns; i++)
                {
                  sBefore = poGetRefln(i)->sGetField(nFI);
                  poGetRefln(i)->vSetField(nFI, sAfter+sBefore);
                }
            }
          else if (1 == nOp)
            {
              for (i = 0; i < m_nNumReflns; i++)
                {
                  sBefore = poGetRefln(i)->sGetField(nFI);
                  poGetRefln(i)->vSetField(nFI, sBefore+sAfter);
                }
            }
          else
            {
              // Illegal operation
              cout << "Creflnlist::nSelect(" << sCompareString << ") is illegal!\n";
              return (-2);
            }
        }
      return (m_nNumReflns);
    }
  else if (eReflnField_int_type == eType)
    {
      return (nSelect(nFI, sOp, nTemp, bIncExc));
    }
  else if (eReflnField_float_type == eType)
    {
      return (nSelect(nFI, sOp, fTemp, bIncExc));
    }
  else if (eReflnField_Cstring_type == eType)
    {
      return (nSelect(nFI, sOp, sAfter, bIncExc));
    }

  return (-3); // Should never get here
}

int
Creflnlist::nExpandRefln(const Creflnlist& oOther, const char* pcSelectIn)
{
  // Add the selected fields in oOther to this list if they do not exist there
  // already.  If pcSelectIn is NULL, then add all fields in the input list.

  int nIntNewFields, nFloatNewFields, nCstringNewFields;
  Cstring *psIntNewNames, *psFloatNewNames, *psCstringNewNames;

  // Loop through fields in oOther and see if they exist in this list.
  // If they do not, check if they are selected, and if so,
  // add them to the items to expand.

  int i;
  int nFI;
  const char *pcSelect;
  const char *pcSelectOther;

  nIntNewFields     = 0;
  nFloatNewFields   = 0;
  nCstringNewFields = 0;
  psIntNewNames     = new Cstring [oOther.m_nIntReflnFields];
  psFloatNewNames   = new Cstring [oOther.m_nFloatReflnFields];
  if (0 < oOther.m_nCstringReflnFields)
    {
      psCstringNewNames = new Cstring [oOther.m_nCstringReflnFields];
    }
  else
    {
      psCstringNewNames = NULL;
    }

  if (NULL == pcSelectIn)
    {
      pcSelect = oOther.m_pcSelect;
    }
  else
    {
      pcSelect = pcSelectIn;
    }

  pcSelectOther = oOther.m_pcSelect;

  // Delete the free list since reflns held there will no longer have matching fields

  (void) nDeleteFree();

  for (i = 0; i < oOther.m_nIntReflnFields; i++)
    {
      psIntNewNames[nIntNewFields] = oOther.sGetFieldName(i, eReflnField_int_type);
      nFI = nGetFieldIndex(psIntNewNames[nIntNewFields], eReflnField_int_type);
      if ( (0 > nFI) && (*pcSelect == *pcSelectOther) )
        {
          // Does not yet exist and not ignored, so add to things to insert

          nIntNewFields++;
        }
      pcSelect++;
      pcSelectOther++;
    }

  for (i = 0; i < oOther.m_nFloatReflnFields; i++)
    {
      psFloatNewNames[nFloatNewFields] = oOther.sGetFieldName(i, eReflnField_float_type);
      nFI = nGetFieldIndex(psFloatNewNames[nFloatNewFields], eReflnField_float_type);
      if ( (0 > nFI) && (*pcSelect == *pcSelectOther) )
        {
          // Does not yet exist and not ignored, so add to things to insert

          nFloatNewFields++;
        }
      pcSelect++;
      pcSelectOther++;
    }

  for (i = 0; i < oOther.m_nCstringReflnFields; i++)
    {
      psCstringNewNames[nCstringNewFields] = oOther.sGetFieldName(i, eReflnField_Cstring_type);
      nFI = nGetFieldIndex(psCstringNewNames[nCstringNewFields], eReflnField_Cstring_type);
      if ( (0 > nFI) && (*pcSelect == *pcSelectOther) )
        {
          // Does not yet exist and not ignored, so add to things to insert

          nCstringNewFields++;
        }
      pcSelect++;
      pcSelectOther++;
    }


  i = nExpandRefln(    nIntNewFields,     psIntNewNames,
                     nFloatNewFields,   psFloatNewNames,
                   nCstringNewFields, psCstringNewNames);

  delete [] psIntNewNames;
  psIntNewNames = NULL;
  delete [] psFloatNewNames;
  psFloatNewNames = NULL;
  if (NULL != psCstringNewNames) 
    {
      delete [] psCstringNewNames;
      psCstringNewNames = NULL;
    }

  return (i);
}

int
Creflnlist::nExpandRefln(const int      nIntNewFields,
                         const Cstring* psIntNewNames,
                         const int      nFloatNewFields,
                         const Cstring* psFloatNewNames,
                         const int      nCstringNewFields,
                         const Cstring* psCstringNewNames)

{
  // Expand the space allocated for reflections in the list by the new number
  // and names of the fields in the argument list

  int i, j;

  int     nIntPrev, nIntNew;
  int     nFloatPrev, nFloatNew;
  int     nCstringPrev, nCstringNew;
  int     nTotalPlus1Prev, nTotalPlus1New;

  Cstring *psPrevNames;
  Crefln  *poRefln;
  char    *pcSelectPrev;

  // Delete the free list since reflns held there will no longer have matching fields

  (void) nDeleteFree();

  // Check if any of the new field names match any old names.
  // If so, that's an error! (May change this in the future!)

  j = 0;
  for (i = 0; i < nIntNewFields; i++)
    {
      if (0 <= nGetFieldIndex(*(psIntNewNames +i), eReflnField_int_type))
        {
          cout << "WARNING expanding reflection list fields: "
               << *(psIntNewNames+i) << " already is a field name!\n";
          j++;
        }
    }
  for (i = 0; i < nFloatNewFields; i++)
    {
      if (0 <= nGetFieldIndex(*(psFloatNewNames +i), eReflnField_float_type))
        {
          cout << "WARNING expanding reflection list fields: "
               << *(psFloatNewNames+i) << " already is a field name!\n";
          j++;
        }
    }
  for (i = 0; i < nCstringNewFields; i++)
    {
      if (0 <= nGetFieldIndex(*(psCstringNewNames +i), eReflnField_Cstring_type))
        {
          cout << "WARNING expanding reflection list fields: "
               << *(psCstringNewNames+i) << " already is a field name!\n";
          j++;
        }
    }
  if (0 < j) return (j);  // Return with error flagged

  // Save number of current reflection fields

  nIntPrev        = m_nIntReflnFields;
  nFloatPrev      = m_nFloatReflnFields;
  nCstringPrev    = m_nCstringReflnFields;
  nTotalPlus1Prev = m_nTotalFieldsPlus1;

  // Set new number of reflection fields

  nIntNew        =     m_nIntReflnFields +     nIntNewFields;
  nFloatNew      =   m_nFloatReflnFields +   nFloatNewFields;
  nCstringNew    = m_nCstringReflnFields + nCstringNewFields;
  nTotalPlus1New = nIntNew + nFloatNew + nCstringNew + 1;

  // Create new pointers to hold all the names (old & new):

  if (0 < nIntNewFields)
    {
      psPrevNames  = m_psIntNames;
      m_psIntNames = new Cstring [ nIntNew ];

      // Copy old names to new list and add the new names

      for (i = 0; i < nIntPrev; i++)             // Copy previous names for Ints
        *(m_psIntNames+i) = *(psPrevNames+i);

      for (j=0; j < nIntNewFields; j++, i++)     // Add new names for Ints
        *(m_psIntNames+i) = *(psIntNewNames+j);

      delete [] psPrevNames;
      psPrevNames = NULL;
    }

  if (0 < nFloatNewFields)
    {
      psPrevNames  = m_psFloatNames;
      m_psFloatNames = new Cstring [ nFloatNew ];

      // Copy old names to new list and add the new names

      for (i = 0; i < nFloatPrev; i++)        // Copy previous names for Floats
        *(m_psFloatNames+i) = *(psPrevNames+i);

      for (j=0; j < nFloatNewFields; j++, i++)   // Add new names for Float
        *(m_psFloatNames+i) = *(psFloatNewNames+j);

      delete [] psPrevNames;
      psPrevNames = NULL;
    }

  if (0 < nCstringNewFields)
    {
      psPrevNames       = m_psCstringNames;
      m_psCstringNames  = new Cstring [ nCstringNew ];

      // Copy old names to new list and add the new names

      for (i = 0; i< nCstringPrev; i++)      // Copy previous names for Cstrings
        *(m_psCstringNames+i) = *(psPrevNames+i);

      for (j=0; j < nCstringNewFields; j++, i++) // Add new names for Cstrings
        *(m_psCstringNames+i) = *(psCstringNewNames+j);

      if (NULL != psPrevNames)
	{
	  delete [] psPrevNames;
	  psPrevNames = NULL;
	}
    }

  pcSelectPrev = m_pcSelect;
  m_pcSelect   = new char [nTotalPlus1New];
  for (i = 0; i< nTotalPlus1New; i++) {
    if (i < nTotalPlus1Prev)
      *(m_pcSelect+i) = *(pcSelectPrev+i);
    else
      *(m_pcSelect+i) = (char) 0;
  }
  delete [] pcSelectPrev;
  pcSelectPrev = NULL;

  // Create a new pointer to point to new list of reflections

  Crefln **ppoPrevList;

  ppoPrevList  = m_ppoTheReflns;
  m_ppoTheReflns = new Crefln* [m_nNumAlloc];

  for (i = 0; i < m_nNumReflns; i++)
    {
      // Allocate new space for reflection, but use new number of fields first

      m_nIntReflnFields     =    nIntNew;
      m_nFloatReflnFields   =  nFloatNew;
      m_nCstringReflnFields = nCstringNew;
      m_nTotalFieldsPlus1   = nTotalPlus1New;

      // In this case, the "free list" cannot be used?!
      // Maybe the "free list" should be deleted by this routine since
      // any reflections there will be the wrong size

      //*(m_ppoTheReflns + i) = new Crefln (this); // Allocate expanded size
      m_ppoTheReflns[i] = new Crefln (this); // Allocate expanded size
      m_nIntReflnFields     =    nIntPrev;       // Restore numfields to prev
      m_nFloatReflnFields   =  nFloatPrev;       //   for the copy coming up
      m_nCstringReflnFields = nCstringPrev;
      m_nTotalFieldsPlus1   = nTotalPlus1Prev;

      poRefln = *(ppoPrevList+i);             // For clarity
      *(*(m_ppoTheReflns + i)) = *poRefln;    // copy Reflection of old size
      delete poRefln;                         // Delete memory for old refln
      poRefln = NULL;
    }

  // Clear old memory holding reflection list

  if (NULL != ppoPrevList)
    {
      delete [] ppoPrevList;
      ppoPrevList = NULL;
    }

  // Reset new number of reflection fields

  m_nIntReflnFields     =    nIntNew;
  m_nFloatReflnFields   =  nFloatNew;
  m_nCstringReflnFields = nCstringNew;
  m_nTotalFieldsPlus1   = nTotalPlus1New;

  vUpdateFieldIndices();                     // Update field indices
  return (0);
}

int
Creflnlist::nGetFieldIndex(const Cstring& sFieldName,
                               const eReflnFieldType eType)
{
  eReflnFieldType eLocalType;

  eLocalType = eType;
  if (eReflnField_unknown_type == eType)
    {
      // If the type is unknown, use the first letter of the sFieldName
      // to try to determine the type.

      if ('n' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_int_type;
        }
      else if ('f' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_float_type;
        }
      else if ('s' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_Cstring_type;
        }
    }

  if (eReflnField_int_type == eLocalType)
    {
      for (int i = 0; i < m_nIntReflnFields; i++)
        if (sFieldName == *(m_psIntNames+i)) return (i);
      return (-1);
    }
  else if (eReflnField_float_type == eLocalType)
    {
      for (int i = 0; i < m_nFloatReflnFields; i++)
        if (sFieldName == *(m_psFloatNames+i)) return (i);
      return (-1);
    }
  else if (eReflnField_Cstring_type == eLocalType)
    {
      for (int i = 0; i < m_nCstringReflnFields; i++)
        if (sFieldName == *(m_psCstringNames+i)) return (i);
      return (-1);
    }
  else
    return(-2);
}

Cstring
Creflnlist::sGetFieldName(const int nFieldIndex,
                          const eReflnFieldType eType) const
{
  if (eReflnField_int_type == eType)
    {
      if ( (nFieldIndex < m_nIntReflnFields) && (nFieldIndex >= 0) )
        return (*(m_psIntNames+nFieldIndex));
    }
  else if (eReflnField_float_type == eType)
    {
      if ( (nFieldIndex < m_nFloatReflnFields) && (nFieldIndex >= 0) )
        return (*(m_psFloatNames+nFieldIndex));
    }
  else if (eReflnField_Cstring_type == eType)
    {
      if ( (nFieldIndex < m_nCstringReflnFields) && (nFieldIndex >= 0) )
        return (*(m_psCstringNames+nFieldIndex));
    }

  return( (Cstring) "");
}

int
Creflnlist::nDelete(const char *pcSelectIn, const bool bIncExc)
{
  // Delete reflections from the reflnlist that are selected or deselected

  const char *pcSelect;
  int        i;                // Loop counter
  int        nNumDeleted = 0;  // Return value (number deleted)

  if (0 >= m_nNumReflns) return (nNumDeleted);  // None to delete;

  if (NULL == pcSelectIn)
    {
      pcSelect = m_pcSelect;
    }
  else
    {
      pcSelect = pcSelectIn;
    }

  int       *pnIndex;

  pnIndex = new int [m_nNumReflns];
  Crefln *poRefln;
  int    nLast = m_nTotalFieldsPlus1 - 1;

  // You know the following if+loops could be easily combined to a single loop

  if (bIncExc)
    {
      // Delete those selected!

      for (i = m_nNumReflns-1; i >=0; i--)
        {
          pnIndex[i] = 0;
          poRefln = *(m_ppoTheReflns+i);
          if (poRefln->m_nSelect == pcSelect[nLast])
            {
              pnIndex[i] = 1;
              nNumDeleted++;
            }
        }
    }
  else
    {
      // Delete those not selected!

      for (i = m_nNumReflns-1; i >=0; i--)
        {
          pnIndex[i] = 0;
          poRefln = *(m_ppoTheReflns+i);
          if (poRefln->m_nSelect != pcSelect[nLast])
            {
              pnIndex[i] = 1;
              nNumDeleted++;
            }
        }
    }

  int nStat = nDelete(1, pnIndex);
  delete [] pnIndex;
  pnIndex = NULL;
  if (nNumDeleted != nStat) 
    cout << "WARNING, inconsistent number of reflections eleted!\n";
  return (nStat);
}

int
Creflnlist::nDelete(const int nFlag, int *pnArray)
{
  // Delete reflns in the array whose value in pnArray[Refnum] == nFlag
  // Returns the number of reflns removed

  int i;
  int *pnTemp;
  int nNumKept = 0;
  int nNumDelt = 0;
  Crefln  *poRefln;
  Crefln **ppoPrevList;
  ppoPrevList  = m_ppoTheReflns;
  m_ppoTheReflns = new Crefln* [m_nNumAlloc];

  pnTemp = pnArray;
  for (i = 0; i < m_nNumReflns; i++, pnTemp++)
    {
      poRefln = ppoPrevList[i];
      if (NULL != poRefln)
        {
          if (nFlag != *pnTemp)
            {
              // Keep it

              m_ppoTheReflns[nNumKept] = poRefln;
              nNumKept++;
            }
          else
            {
              nNumDelt++;
              // Instead of actually deleting the poRefln,
              // it is better to put it on a "free list"

              if (m_nNumAllocFree <= m_nNumFree)
                {
                  // Create space in the free list for this "deleted" reflection
        
                  int      j;
                  Crefln **ppoPrevFreeList;
                  ppoPrevFreeList = m_ppoFreeReflns;
                  m_nNumAllocFree = m_nNumAllocFree + m_nAllocSize;
                  m_ppoFreeReflns = new Crefln* [m_nNumAllocFree];

                  for (j = 0; j < m_nNumFree; j++)
                    {
                      m_ppoFreeReflns[j] = ppoPrevFreeList[j];
                    }
                  for (j = m_nNumFree; j < m_nNumAllocFree; j++)
                    m_ppoFreeReflns[j] = NULL;
                  delete [] ppoPrevFreeList;
		  ppoPrevFreeList = NULL;
                }

              // Place this "deleted" reflection in the free list

              m_ppoFreeReflns[m_nNumFree] = poRefln;
              m_nNumFree++;
            }
        }
      else
        {
          cout << "WARNING! nDelete NULL refln i= " << i << '\n' << flush;
        }
    }

  // Some other clean-up

  if ( (NULL != m_pnSortIndex) && (0 < nNumDelt) )
    {
      delete [] m_pnSortIndex;
      m_pnSortIndex = NULL;
    }

  int nStat;
  nStat = m_nNumReflns - nNumKept;
  m_nNumReflns = nNumKept;
  for (i = m_nNumReflns; i < m_nNumAlloc; i++)
    m_ppoTheReflns[i] = NULL;

  if (m_nNumReflns==0)
      nEraseExtra();

  delete [] ppoPrevList;
  ppoPrevList = NULL;
  return (nStat);
}

int
Creflnlist::nDelete(const int nReflnNum) // Delete a reflection
{
  Crefln *poRefln;

  if ( (0 <= nReflnNum) && (nReflnNum < m_nNumReflns) )
    {
      poRefln = *(m_ppoTheReflns+nReflnNum);
      if (NULL != poRefln)
        {
          // Instead of actually deleting the poRefln,
          // it would be better to put it on a "free list"

//          delete poRefln;

          if (m_nNumAllocFree <= m_nNumFree)
            {
              // Create space in the free list for this "deleted" reflection
        
              int      i;
              Crefln **ppoPrevList;
              ppoPrevList  = m_ppoFreeReflns;
              m_nNumAllocFree = m_nNumAllocFree + m_nAllocSize;
              m_ppoFreeReflns = new Crefln* [m_nNumAllocFree];

              for (i = 0; i < m_nNumFree; i++)
                {
                  m_ppoFreeReflns[i] = ppoPrevList[i];
                }
              for (i = m_nNumFree; i < m_nNumAllocFree; i++)
                m_ppoFreeReflns[i] = NULL;
              if (ppoPrevList)
		{
                  delete [] ppoPrevList;
		  ppoPrevList = NULL;
		}
            }

          // Place this "deleted" reflection in the free list

          m_ppoFreeReflns[m_nNumFree] = poRefln;
          m_nNumFree++;
        }
      m_nNumReflns--;
      for (int i = nReflnNum; i < m_nNumReflns; i++)
        {
          *(m_ppoTheReflns+i) = *(m_ppoTheReflns+i+1);// Shift all ptrs down one
        }
      *(m_ppoTheReflns+m_nNumReflns) = NULL;          // Set last to 0.

      if (NULL != m_pnSortIndex)
        {
          delete [] m_pnSortIndex;
          m_pnSortIndex = NULL;
        }

      if (m_nNumReflns==0)
        nEraseExtra();

      return (0);
    }
  else
    {
      cout << "Reflection index " << nReflnNum << " not in reflection list!\n";
      return (1);
    }
}

/* nDeleteSorted() deletes a reflection, but also adjusts the m_pnSortIndex[] array.
   Also, the user can provide various external arrays that must also be adjusted for the deletion.

  */

int Creflnlist::nDeleteSorted(const int nReflnNumSorted1,
                              int* pnArray1, 
                              int* pnArray2, 
                              int* pnArray3) 
{
    
    int* pnSortIndex;
    int nSortArray;
    int* pnIndex;
    int i;
    int nReflnNumReal;
    int nReflnNumSorted;

    if (!pnArray1)
        return 1;
    nReflnNumReal = pnArray1[nReflnNumSorted1];
    
    for (nSortArray=1;nSortArray<4;nSortArray++) {

        if (nSortArray==0) 
            pnIndex = m_pnSortIndex;
        else if (nSortArray==1)
            pnIndex = pnArray1;
        else if (nSortArray==2)
            pnIndex = pnArray2;
        else if (nSortArray==3)
            pnIndex = pnArray3;
        else
            pnIndex = NULL;

        if (pnIndex) {

            if (nSortArray==1)
                nReflnNumSorted = nReflnNumSorted1;
            else {
                // Must search for the 'real' index w.r.t. this sort array.
                for (i=0;i< nGetNumReflns(); i++) {
                    if (pnIndex[i] == nReflnNumReal)
                        break;
                }
                if (i == nGetNumReflns())
                    return 1;
                nReflnNumSorted = i;
            }

            for (i = 0; i< nGetNumReflns() - 1 ; i++) {
                if (i>=nReflnNumSorted) {
                    pnIndex[i] = pnIndex[i + 1];
                }
                if (pnIndex[i]>nReflnNumReal)
                    pnIndex[i]--;
            }
        }
    }

    pnSortIndex = m_pnSortIndex;
    m_pnSortIndex = NULL;       // Save the sort index.  Otherwise, nDelete() deletes it!
    nDelete(nReflnNumReal);
    m_pnSortIndex = pnSortIndex;

    return 0;
}

void
Creflnlist::vDeleteAll(void)
{
  // Delete all the reflections in the list
  // (note they may go into the

  int i;
  int nNumReflns;
  nNumReflns = nGetNumReflns();
  for (i = nNumReflns-1; i >=0; i--)
    nDelete(i);
  if (m_nNumReflns==0)
      nEraseExtra();
}

int
Creflnlist::nExpandGetField(const Cstring&  sFieldName,
                                const eReflnFieldType eType)
{
  // Return the field index of a field.  If it doesn't exist, create the
  // field and expand the reflection list.  Return value > 0 on error

  int nIndex;
  eReflnFieldType eLocalType;

  eLocalType = eType;
  if (eReflnField_unknown_type == eLocalType)
    {
      // If the type is unknown, use the first letter of the sFieldName
      // to try to determine the type.

      if ('n' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_int_type;
        }
      else if ('f' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_float_type;
        }
      else if ('s' == sFieldName.GetAt(0))
        {
          eLocalType = eReflnField_Cstring_type;
        }
    }

  // Now look for the index

  nIndex = nGetFieldIndex(sFieldName, eLocalType);

  if (-1 == nIndex)
    {
      // Field does not exist in the reflection list, so try to add it

      if (eReflnField_int_type == eLocalType)
        {
          nIndex = nExpandRefln(1, &sFieldName);
        }
      else if (eReflnField_float_type == eLocalType)
        {
          nIndex = nExpandRefln(0, NULL, 1, &sFieldName);
        }
      else if (eReflnField_Cstring_type == eLocalType)
        {
          nIndex = nExpandRefln(0, NULL, 0, NULL, 1, &sFieldName);
        }
      if (0 == nIndex)
        // Check return status of nExpandRefln(...), if success, get the index
        {
          nIndex = nGetFieldIndex(sFieldName, eLocalType);
        }
      else if (0 < nIndex)
        {
          nIndex = -nIndex;  // Convert error to < 0
        }
    }
  vUpdateFieldIndices();   // Update any member field indices
  return (nIndex);
}

int
Creflnlist::nInsertReflnFrom(Crefln& oReflnNative,Crefln& oReflnOther) {
    int i;

    // Ok, this reflection is NOT excluded, so it can be inserted
    for (i = 0; i < m_nIntReflnFields; i++)
    {
        if (0 <= m_pnIntFieldsRelative[i])
            oReflnNative.vSetField(i, oReflnOther.nGetField(m_pnIntFieldsRelative[i]));
    }
    
    for (i = 0; i < m_nFloatReflnFields; i++)
    {
        if (0 <= m_pnFloatFieldsRelative[i])
            oReflnNative.vSetField(i, oReflnOther.fGetField(m_pnFloatFieldsRelative[i]));
    }
    
    for (i = 0; i < m_nCstringReflnFields; i++)
    {
        if (0 <= m_pnCstringFieldsRelative[i])
            oReflnNative.vSetField(i, oReflnOther.sGetField(m_pnCstringFieldsRelative[i]));
    }
    
    // OK, oRefln has every matching item transferred, so insert into list.
    
    return nInsert(&oReflnNative);
}


int
Creflnlist::nInsertListFrom(Creflnlist& oReflnlist, const char *pcSelectIn, const int* pnSortOrder,int nNumToInsert)
{
  // Insert all reflections from another list, but only the fields that
  //   are already in this list.  Do not add any new fields to this list.

  // ToDo:  Maybe insert only the reflections that are selected!

  int     i, j;             // Loop counter
  int     nStat;            // General status value
  Crefln  oRefln(this);     // A reflection object for filling and inserting
  Crefln *poOtherRefln;     // A reflection from the other list.

  const char *pcSelect;           // Selection field

//  int nIntFields[nIntReflnFields];
//  int nFloatFields[nFloatReflnFields];
//  int nCstringFields[nCstringReflnFields];

  if (m_pnIntFieldsRelative)
      delete[] m_pnIntFieldsRelative;
  m_pnIntFieldsRelative     = new int [m_nIntReflnFields];
  if (m_pnFloatFieldsRelative)
      delete[] m_pnFloatFieldsRelative;
  m_pnFloatFieldsRelative     = new int [m_nFloatReflnFields];
  if (m_pnCstringFieldsRelative)
      delete[] m_pnCstringFieldsRelative;
  m_pnCstringFieldsRelative     = new int [m_nCstringReflnFields];

  // Setup selection mechanism.  Only reflections whose
  // selection field matches the input will be inserted into the new list.

  if (NULL == pcSelectIn)
    {
      pcSelect = oReflnlist.m_pcSelect;
    }
  else
    {
      pcSelect = pcSelectIn;
    }

  if (nNumToInsert == -1)
      nNumToInsert = oReflnlist.nGetNumReflns();

  // To make things faster, make sure enough slots are available for
  // the new reflns.

  nExpand(nNumToInsert + 200);

  // Get matching field names and initialize matching field indices
  // If the field name from this list exists in the input list,
  // then the index will be >= 0.  If it does not exist, then the index
  // will be < 0.

  for (i = 0; i < m_nIntReflnFields; i++)
    {
      m_pnIntFieldsRelative[i] = oReflnlist.nGetFieldIndex(sGetFieldName(i,
                                                    eReflnField_int_type),
                                                eReflnField_int_type);
    }

  for (i = 0; i < m_nFloatReflnFields; i++)
    {
      m_pnFloatFieldsRelative[i] = oReflnlist.nGetFieldIndex(sGetFieldName(i,
                                                      eReflnField_float_type),
                                                  eReflnField_float_type);
    }

  for (i = 0; i < m_nCstringReflnFields; i++)
    {
      m_pnCstringFieldsRelative[i] = oReflnlist.nGetFieldIndex(sGetFieldName(i,
                                                        eReflnField_Cstring_type),
                                                    eReflnField_Cstring_type);
    }

  // Loop through reflections, copy matching fields, then insert


  nStat = 0;
  for (j = 0; ( (j < nNumToInsert) && (0 == nStat) ); j++)
    {
      if (pnSortOrder)
          poOtherRefln = oReflnlist.m_ppoTheReflns[pnSortOrder[j]];
      else
          poOtherRefln = oReflnlist.m_ppoTheReflns[j];

      if (   pcSelect[oReflnlist.m_nTotalFieldsPlus1-1]
          == poOtherRefln->m_nSelect)
        {
          nStat = nInsertReflnFrom(oRefln,*poOtherRefln);
        }
    }

  vUpdateFieldIndices();
  nAddExtra(oReflnlist);

  return (nStat);
}

void 
Creflnlist::vSort2(const eReflnFieldType eType1, const int nField1,
                   const eReflnFieldType eType2, const int nField2,
                   int* pnIndexIn) 
{
  int* pnIndex;
  int* pnBufPt;
  int  nx,ny,nz;

  if (NULL != pnIndexIn) 
    {
     pnIndex = pnIndexIn;
   } 
  else 
    {
      if (NULL != m_pnSortIndex)
        delete [] m_pnSortIndex;
      m_pnSortIndex = new int [max(1,nGetNumReflns()+1)];
      m_pnSortIndex[nGetNumReflns()] = 0;
      pnIndex       = m_pnSortIndex;
    }
  pnBufPt = new int[max(1,nGetNumReflns())];

  if (0 < nGetNumReflns())
    {
     if (eReflnField_float_type == eType1) 
       {
        vSortFloatHeap(nField1, pnIndex);

        // Now, search for equivalent reflections, and sort sub lists.
        for (nx=1, ny=0; nx <= nGetNumReflns(); nx++)
          {
            if  (   (nx == nGetNumReflns()) 
                 || (!(m_ppoTheReflns[pnIndex[nx]]->fGetField(nField1)
                       == m_ppoTheReflns[pnIndex[nx-1]]->fGetField(nField1)))) 
              {
                if (nx-ny <= 1) {
                  ny=nx;
                  continue;
                }
                for (nz=0; nz < nx-ny; nz++)
                  pnBufPt[nz]=pnIndex[ny+nz];
                if (eReflnField_float_type == eType2)
                  vSortFloatHeap(nField2, pnBufPt, 0, nx-ny, TRUE);
                else if (eReflnField_int_type == eType2)
                  vSortIntHeap(nField2, pnBufPt, 0, nx-ny, TRUE);
                else if (eReflnField_Cstring_type == eType2)
                  vSortStringRadixRightAdjusted(nField2, pnBufPt, 0, nx-ny, TRUE);
                for (nz=0; nz < nx-ny; nz++) 
                  pnIndex[ny+nz] = pnBufPt[nz];
                ny = nx;
              }
          }
      } 
     else if (eReflnField_int_type == eType1)
       {
         vSortIntHeap(nField1, pnIndex);
         
         // Now, search for equivalent reflections, and sort sub lists.
         
         for (nx=1,ny=0; nx <= nGetNumReflns(); nx++) 
           {
             if (   (nx == nGetNumReflns()) 
                 || (!(m_ppoTheReflns[pnIndex[nx]]->nGetField(nField1)
                       == m_ppoTheReflns[pnIndex[nx-1]]->nGetField(nField1)))) 
               {
                 if (nx-ny <= 1) {
                   ny=nx;
                   continue;
                 }
                 for (nz = 0;nz < nx-ny; nz++) 
                   pnBufPt[nz] = pnIndex[ny+nz];
                 if (eReflnField_float_type == eType2)
                   vSortFloatHeap(nField2,pnBufPt, 0, nx-ny, TRUE);
                 else if (eReflnField_int_type == eType2)
                   vSortIntHeap(nField2, pnBufPt, 0, nx-ny, TRUE);
                 else if (eReflnField_Cstring_type == eType2) 
                   vSortStringRadixRightAdjusted(nField2, pnBufPt, 0, nx-ny, TRUE);
                 for (nz=0; nz < nx-ny; nz++) 
                   pnIndex[ny+nz] = pnBufPt[nz];
                 ny = nx;
               }
           }  
       } 
     else if (eReflnField_Cstring_type == eType1) 
       {
         vSortStringRadixRightAdjusted(nField1, pnIndex);

        // Now, search for equivalent reflections, and sort sub lists.
        for (nx=1,ny=0;nx<=nGetNumReflns();nx++) 
          {
            if  (   (nx == nGetNumReflns())
                 || (!(m_ppoTheReflns[pnIndex[nx]]->sGetField(nField1)
                       == m_ppoTheReflns[pnIndex[nx-1]]->sGetField(nField1)))) 
              {
               if (nx-ny <= 1) {
                    ny=nx;
                    continue;
               }
               for (nz=0; nz < nx-ny; nz++) 
                 pnBufPt[nz] = pnIndex[ny+nz];
               if (eReflnField_float_type == eType2)
                 vSortFloatHeap(nField2, pnBufPt, 0, nx-ny, TRUE);
               else if (eReflnField_int_type == eType2)
                 vSortIntHeap(nField2, pnBufPt, 0, nx-ny, TRUE);
               else if (eReflnField_Cstring_type == eType2)
                 vSortStringRadixRightAdjusted(nField2, pnBufPt, 0, nx-ny, TRUE);
               for (nz=0; nz < nx-ny; nz++) 
                 pnIndex[ny+nz] = pnBufPt[nz];
               ny = nx;
             }
          }
      }
   }

  delete[] pnBufPt;
  return;
}

////////////////////////////////////////////////////////////
void Creflnlist::vSort(const eReflnFieldType eType, 
                       const int nField,
                       int* pnIndexIn)
{
    int*        pnIndex = NULL;
    
    if( NULL != pnIndexIn )
    {
        pnIndex = pnIndexIn;
    }
    else
    {
        if( NULL != m_pnSortIndex )
            delete [] m_pnSortIndex;

        m_pnSortIndex = new int [max(1,nGetNumReflns()+1)];
        
        m_pnSortIndex[nGetNumReflns()] = 0;
        
        pnIndex = m_pnSortIndex;
    }
    
    if( 0 < nGetNumReflns() )
    {
        if( eReflnField_float_type == eType )
            vSortFloatHeap(nField, pnIndex);
        else if( eReflnField_int_type == eType )
            vSortIntHeap(nField, pnIndex);
        else if( eReflnField_Cstring_type == eType )
            vSortStringRadixRightAdjusted(nField,pnIndex);
    }
}
////////////////////////////////////////////////////////////////
void
Creflnlist::vSortFloatHeap(const int nField, int *pnIndex,int nFirstRefln,int nNumReflns,bool bUseIndexOrder)
{
  // Heap sort on index from Numerical Recipes... Section 8.3

  //  cout << "Sorting on field: " << m_psIntNames[nField] << endl;
  int   indxt, ir, i, j, m;
  float fTemp;
  int  *pnTemp;

  if (0 == nNumReflns) 
    {
      if (bUseIndexOrder) 
        return;
      nNumReflns = m_nNumReflns-nFirstRefln;
    }

  if ((bUseIndexOrder) && (0 != nFirstRefln)) 
    {
      cout << "Mode Not supported yet in vSort????Heap()\n" << flush;
      exit(1);   // Programmer error, watch out SSI folks!
    }

  if (!bUseIndexOrder) 
    {
      for (i = 0; i < nNumReflns; i++) 
        pnIndex[i] = i+nFirstRefln;  // Init w/consec integers
    }

  if (1 >= nNumReflns)
    {
      return;                      // That was easy!
    }


  pnTemp = pnIndex-1;        // A little kludge per NR for C arrays, so
                             // indexing starts from 1
  m  = nNumReflns / 2 + 1;
  ir = nNumReflns;

  for (;;)
    {
      if (1 < m)
        {
          m--;
          indxt = pnTemp[m];
          fTemp = m_ppoTheReflns[indxt]->fGetField(nField);
        }
      else
        {
          indxt       = pnTemp[ir];
          fTemp       = m_ppoTheReflns[indxt]->fGetField(nField);
          pnTemp[ir] = pnTemp[1];
          ir--;
          if (1 == ir)
            {
              pnTemp[1] = indxt;
              break;
              ///////!!!!
            }
        }
      i = m;
      j = m + m;
      while (j <= ir)
        {
          if (   (j < ir)
              && (m_ppoTheReflns[pnTemp[j]]->fGetField(nField)
                < m_ppoTheReflns[pnTemp[j+1]]->fGetField(nField)) )
            j++;

          if (fTemp < m_ppoTheReflns[pnTemp[j]]->fGetField(nField))
            {
              pnTemp[i] = pnTemp[j];
              i = j;
              j = j + j;
            }
          else
            {
              j = ir + 1;
            }
        }
      pnTemp[i] = indxt;
    } // end of for ;;
}

void
Creflnlist::vSortIntRadix(const int nField, int *pnIndex,
                          int nFirstRefln, int nNumReflns) {
    const int nRadixShift = 2;
    const int nRadixArrays = 1 << nRadixShift;
    const unsigned int unAnd = nRadixArrays - 1;
    const int nIntBits = sizeof(int)*8;

    int* apnIndex[nRadixArrays];
    int  anIndexUsed[nRadixArrays];
    int nx,ny;
    unsigned int unMaxInt = 1;
    int nShift;

    if (0 == nNumReflns) {
        nNumReflns = m_nNumReflns - nFirstRefln;
    }

    for (nx = 0; nx < nNumReflns; nx++) 
        pnIndex[nx] = nx+nFirstRefln;  // Init w/consec integers

    if (1 >= nNumReflns) 
    {
      return;                      // That was easy!
    }

    // Allocate radix arrays.
    for (nx=0;nx<nRadixArrays;nx++) {
        apnIndex[nx] = new int[nNumReflns];
    }

    // We make passes until we have gone over all 32 bits of an "int" or we have gone past the largest integer.
    for (nShift=0; (nShift < nIntBits) && ((((unsigned int) 1) << nShift)<=unMaxInt); nShift += nRadixShift) {
        
        // Reset the maximum integer.
        unMaxInt=0;

        // Set all of the arrays so that they have no data in them.
        for (nx=0;nx<nRadixArrays;nx++) 
            anIndexUsed[nx] = 0;

        int nRef;
        int nIndex;
        int nArray;
        unsigned int unValue;

        // Sort the values into the arrays.

        for (nRef = 0; nRef<nNumReflns; nRef++) {
            // It is a little tricky casting to an 
            unValue = (unsigned int) m_ppoTheReflns[pnIndex[nRef]]->nGetField(nField);
            nArray = unAnd & (unValue >> nShift);
            apnIndex[nArray][anIndexUsed[nArray]++] = pnIndex[nRef];           
            if (unValue>unMaxInt)
                unMaxInt = unValue;
        }

        // Replace the values into pnIndex array.
        
        nIndex = 0;
        for (nx = 0; nx<nRadixArrays; nx++) {
            for (ny=0;ny<anIndexUsed[nx]; ny++)
                pnIndex[nIndex++] = apnIndex[nx][ny];           
        }
    }


    // De-Allocate radix arrays.
    for (nx=0;nx<nRadixArrays;nx++) {
        delete[] apnIndex[nx];
    }
    return;
}


void Creflnlist::vSortStringRadix(const int nField, 
                                  int *pnIndex,
                                  int nFirstRefln, 
                                  int nNumReflns) 
{

    const int nRadixArrays = nCREFLN_STRING_SIZE-1;

    int         nx = 0;
    int         ny = 0;
    int         nz = 0;

    int         nRef = 0;
    char*       pcTemp = NULL;

    bool        abCharUsed[256];


    if( 0 == nNumReflns )
    {
        nNumReflns = m_nNumReflns - nFirstRefln;
    }

    for(nx=0; nx < nNumReflns; nx++)
        pnIndex[nx] = nx + nFirstRefln;

    if( 1 >= nNumReflns ) 
    {
        return;                      // That was easy!
    }

    int         nLength = 0;
    int         nLengthBound = 0;

    for(nRef = 0; nRef < nNumReflns; nRef++)
    {
        pcTemp = (char*)m_ppoTheReflns[pnIndex[nRef]]->pcGetField(nField);
        
        nLength = strlen(pcTemp);
        
        nLengthBound = max(nLength,nLengthBound);
    }

    // Allocate radix arrays.
    int*        pnIndex256 = new int [nNumReflns];
    int*        pnIndexNext = new int [nNumReflns];
    
    // We make passes until we have gone over all 32 bits of an "int" or we have gone past the largest integer.
    for(int nShift = nLengthBound - 1; nShift >= 0; nShift--)
    {
        // Sort the values into the arrays.
        memset(abCharUsed,0,sizeof(bool)*256);

        for (nRef = 0; nRef < nNumReflns; nRef++)
        {
            pcTemp = (char*)m_ppoTheReflns[pnIndex[nRef]]->pcGetField(nField);
            nLength = strlen(pcTemp);
            
            if( nLength <= nShift)
                pnIndex256[nRef] = 0;
            else
                pnIndex256[nRef] = (int)pcTemp[nShift];
            
            abCharUsed[pnIndex256[nRef]] = 1;            
        }
        
        for(nx = 0, nz = 0; nx < 256; nx++)
        {
            if( abCharUsed[nx] )
            {
                for(ny = 0; ny < nNumReflns; ny++)
                {
                    if( pnIndex256[ny] == nx )
                    {
                        pnIndexNext[nz++] = pnIndex[ny];
                    }
                }
            }
        }
        
        for(nx = 0; nx < nNumReflns; nx++)
            pnIndex[nx] = pnIndexNext[nx];
    }

    if( pnIndex256 )
        delete[] pnIndex256;
    
    if( pnIndexNext )
        delete[] pnIndexNext;
}
////////////////////////////////////////////////////////////////////
void Creflnlist::vSortStringRadixRightAdjusted(const int nField, 
                                               int *pnIndex,
                                               int nFirstRefln, 
                                               int nNumReflns,
                                               bool bUseInputIndexArray) 
{
    const int nRadixArrays = nCREFLN_STRING_SIZE-1;

    int         nx = 0;
    int         ny = 0;
    int         nz = 0;

    int         nRef = 0;
    char*       pcTemp = NULL;

    bool        abCharUsed[256];

    if( 0 == nNumReflns )
    {
        nNumReflns = m_nNumReflns - nFirstRefln;
    }

    if( !bUseInputIndexArray )
    {
        for(nx=0; nx < nNumReflns; nx++)
            pnIndex[nx] = nx + nFirstRefln;
    }

    if( 1 >= nNumReflns ) 
    {
        return;                      // That was easy!
    }

    int         nLength = 0;
    int         nLengthBound = 0;

    for(nRef = 0; nRef < nNumReflns; nRef++)
    {
        pcTemp = (char*)m_ppoTheReflns[pnIndex[nRef]]->pcGetField(nField);
        
        nLength = strlen(pcTemp);
        
        nLengthBound = max(nLength,nLengthBound);
    }

    // Allocate radix arrays.
    int*        pnIndex256 = new int [nNumReflns];
    int*        pnIndexNext = new int [nNumReflns];
    
    // We make passes until we have gone over all 32 bits of an "int" or we have gone past the largest integer.
    for(int nShiftFromEndToLeft = 0; nShiftFromEndToLeft < nLengthBound; nShiftFromEndToLeft++)
    {
        // Sort the values into the arrays.
        memset(abCharUsed, false, sizeof(bool)*256);

        for (nRef = 0; nRef < nNumReflns; nRef++)
        {
            pcTemp = (char*)m_ppoTheReflns[pnIndex[nRef]]->pcGetField(nField);
            nLength = strlen(pcTemp);
            
            if( nLength <= nShiftFromEndToLeft)
                pnIndex256[nRef] = 0;
            else
                pnIndex256[nRef] = (int)pcTemp[nLength - nShiftFromEndToLeft - 1];
            
            abCharUsed[pnIndex256[nRef]] = 1;            
        }
        
        for(nx = 0, nz = 0; nx < 256; nx++)
        {
            if( abCharUsed[nx] )
            {
                for(ny = 0; ny < nNumReflns; ny++)
                {
                    if( pnIndex256[ny] == nx )
                    {
                        pnIndexNext[nz++] = pnIndex[ny];
                    }
                }
            }
        }
        
        for(nx = 0; nx < nNumReflns; nx++)
            pnIndex[nx] = pnIndexNext[nx];
    }

    if( pnIndex256 )
        delete[] pnIndex256;
    
    if( pnIndexNext )
        delete[] pnIndexNext;
}



void
Creflnlist::vSortStringHeap(const int nField, int *pnIndex,
                            int nFirstRefln, int nNumReflns, 
                            bool bUseIndexOrder)
{
  // Heap sort on index from Numerical Recipes... Section 8.3

  //  cout << "Sorting on field: " << m_psIntNames[nField] << endl;
  int   indxt, ir, i, j, m;
  Cstring sTemp;
  int  *pnTemp;

  if (0 == nNumReflns)
    {
      if (bUseIndexOrder) 
        return;
      nNumReflns = m_nNumReflns-nFirstRefln;
  }

  if ((bUseIndexOrder) && (0 != nFirstRefln))
    { 
      cout << "Mode Not supported yet in vSort????Heap()\n" << flush;
      exit(1);
    }

  if (!bUseIndexOrder) 
    {
      for (i = 0; i < nNumReflns; i++) 
        pnIndex[i] = i+nFirstRefln;  // Init w/consec integers
    }

  if (1 >= nNumReflns) 
    {
      return;                      // That was easy!
    }

  pnTemp = pnIndex-1;        // A little kludge per NR for C arrays, so
                             // indexing starts from 1
  m  = nNumReflns / 2 + 1;
  ir = nNumReflns;

  for (;;)
    {
      if (1 < m)
        {
          m--;
          indxt = pnTemp[m];
          sTemp = m_ppoTheReflns[indxt]->sGetField(nField);
        }
      else
        {
          indxt       = pnTemp[ir];
          sTemp       = m_ppoTheReflns[indxt]->sGetField(nField);
          pnTemp[ir] = pnTemp[1];
          ir--;
          if (1 == ir)
            {
              pnTemp[1] = indxt;
              break;
              ///////!!!!
            }
        }
      i = m;
      j = m + m;
      while (j <= ir)
        {
          if (   (j < ir)
              && (m_ppoTheReflns[pnTemp[j]]->sGetField(nField)
                < m_ppoTheReflns[pnTemp[j+1]]->sGetField(nField)) )
            j++;

          if (sTemp < m_ppoTheReflns[pnTemp[j]]->sGetField(nField))
            {
              pnTemp[i] = pnTemp[j];
              i = j;
              j = j + j;
            }
          else
            {
              j = ir + 1;
            }
        }
      pnTemp[i] = indxt;
    } // end of for ;;

}


void
Creflnlist::vSortIntHeap(const int nField, int* pnIndex,int nFirstRefln,int nNumReflns,bool bUseIndexOrder)
{
  // Heap sort on index from Numerical Recipes... Section 8.3

  //  cout << "Sorting on field: " << m_psIntNames[nField] << endl;
  int   indxt, ir, i, j, m;
  int   nTemp;
  int  *pnTemp;

  if (0 == nNumReflns)
    {
      if (bUseIndexOrder) 
        return;
      nNumReflns = m_nNumReflns-nFirstRefln;
  }

  if ((bUseIndexOrder) && (0 != nFirstRefln))
    {
      cout << "Mode Not supported yet in vSort????Heap()\n" << flush;
      exit(1);
    }

  if (!bUseIndexOrder) 
    {
      for (i = 0; i < nNumReflns; i++) 
        pnIndex[i] = i+nFirstRefln;  // Init w/consec integers
    }

  if (1 >= nNumReflns)
    {
      return;                      // That was easy!
    }

  pnTemp = pnIndex-1;        // A little kludge per NR for C arrays, so
                             // indexing starts from 1
  m  = nNumReflns / 2 + 1;
  ir = nNumReflns;

  for (;;)
    {
      if (1 < m)
        {
          m--;
          indxt = pnTemp[m];
          nTemp = m_ppoTheReflns[indxt]->nGetField(nField);
        }
      else
        {
          indxt       = pnTemp[ir];
          nTemp       = m_ppoTheReflns[indxt]->nGetField(nField);
          pnTemp[ir] = pnTemp[1];
          ir--;
          if (1 == ir)
            {
              pnTemp[1] = indxt;
              break;
              ///////!!!!
            }
        }
      i = m;
      j = m + m;
      while (j <= ir)
        {
          if (   (j < ir)
              && (m_ppoTheReflns[pnTemp[j]]->nGetField(nField)
                < m_ppoTheReflns[pnTemp[j+1]]->nGetField(nField)) )
            j++;

          if (nTemp < m_ppoTheReflns[pnTemp[j]]->nGetField(nField))
            {
              pnTemp[i] = pnTemp[j];
              i = j;
              j = j + j;
            }
          else
            {
              j = ir + 1;
            }
        }
      pnTemp[i] = indxt;
    } // end of for ;;

}

void
Creflnlist::vSortAbsFloatHeap(const int nField, int *pnIndex)
{
  // Heap sort on index from Numerical Recipes... Section 8.3

  //  cout << "Sorting on field: " << m_psIntNames[nField] << endl;

  int   indxt, ir, i, j, m;
  float fTemp;
  int  *pnTemp;

  for (i = 0; i < m_nNumReflns; i++) 
    pnIndex[i] = i;  // Init w/consec integers

  if (1 >= m_nNumReflns) 
    {
      return;                      // That was easy!
    }

  pnTemp = pnIndex-1;        // A little kludge per NR for C arrays, so
                             // indexing starts from 1
  m  = m_nNumReflns / 2 + 1;
  ir = m_nNumReflns;

  for (;;)
    {
      if (1 < m)
        {
          m--;
          indxt = pnTemp[m];
          fTemp = m_ppoTheReflns[indxt]->fGetField(nField);
          fTemp = (float)fabs((double)fTemp);
        }
      else
        {
          indxt       = pnTemp[ir];
          fTemp       = m_ppoTheReflns[indxt]->fGetField(nField);
          fTemp       = (float)fabs((double)fTemp);
          pnTemp[ir] = pnTemp[1];
          ir--;
          if (1 == ir)
            {
              pnTemp[1] = indxt;
              return;
              ///////!!!!
            }
        }
      i = m;
      j = m + m;
      while (j <= ir)
        {
          if (   (j < ir)
              && (fabs((double)m_ppoTheReflns[pnTemp[j]]->fGetField(nField))
                < fabs((double)m_ppoTheReflns[pnTemp[j+1]]->fGetField(nField))) )
            j++;

          if (fTemp < (float)fabs((double)m_ppoTheReflns[pnTemp[j]]->fGetField(nField)))
            {
              pnTemp[i] = pnTemp[j];
              i = j;
              j = j + j;
            }
          else
            {
              j = ir + 1;
            }
        }
      pnTemp[i] = indxt;
    } // end of for ;;
}


/*  These routines execute a binary search for values.
*/

int  Creflnlist::nFindFirst(int nField,const int nValue,int* pnIndex,bool bFindFirstGreater) {
    int* pnSortIndex;
    int  nFirst,nLast,nTest;
    int  nFirstValue,nLastValue,nTestValue;
    int  nPass;

    if (nGetNumReflns()==0) {
        if (bFindFirstGreater)
            return 0;
        else
            return -1;
    }

    if (pnIndex)
        pnSortIndex = pnIndex;
    else
        pnSortIndex = m_pnSortIndex;

    nFirst = 0;
    nLast = nGetNumReflns()-1;
    nPass = 0;

    while (1) {
        nFirstValue = (*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].nGetField(nField);
        nLastValue = (*this)[(pnSortIndex)?(pnSortIndex[nLast]):nLast].nGetField(nField);

        if ((nPass==0) && (nLastValue<nValue)) {
            // If we are finding the first greater, then this reflection gets put at the end of the list.
            if (bFindFirstGreater)
                return m_nNumReflns;
            else
                return -1;
        }

        if ((nPass==0) && (nFirstValue>nValue)) {
            // If we are finding the first greater, then this reflection gets put at the beginning of the list.
            if (bFindFirstGreater)
                return 0;
            else
                return -1;
        }
        if (nFirst==nLast)
            break;
        
        if (nFirstValue==nValue)
            break;

        nTest = (nFirst+nLast)/2+1;

        nTestValue = (*this)[(pnSortIndex)?(pnSortIndex[nTest]):nTest].nGetField(nField);
        if (nTestValue>nValue)
            nLast = nTest-1;
        else 
            nFirst = nTest;

        nPass++;
    }
    if (bFindFirstGreater) {
        // We know have *[nFirst] <= nValue 
        // We need to look forward to the first element that is > nValue.
        while ((nFirst+1<m_nNumReflns) && ((*this)[(pnSortIndex)?(pnSortIndex[nFirst+1]):(nFirst+1)].nGetField(nField)<=nValue))
            nFirst++;
        if (nFirst<m_nNumReflns)
            nFirst++;
        return nFirst;
    } else {
        // did we find the value?    
        if ((nFirst>=nGetNumReflns()) || ((*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].nGetField(nField)!=nValue)) 
            return -1;
        // Yes we did...
        // Now backtrack to find the first element that is == nValue
        while ((nFirst>0) && ((*this)[(pnSortIndex)?(pnSortIndex[nFirst-1]):(nFirst-1)].nGetField(nField)==nValue))
            nFirst--;
        return nFirst;
    }
}

int  Creflnlist::nFindFirst(int nField,const double fValue,double fTol,int* pnIndex,bool bFindFirstGreater) {
    int* pnSortIndex;
    int  nFirst,nLast,nTest;
    double  fFirstValue,fLastValue,fTestValue;
    int  nPass;

    if (nGetNumReflns()==0) {
        if (bFindFirstGreater)
            return 0;
        else
            return -1;
    }

    if (pnIndex)
        pnSortIndex = pnIndex;
    else
        pnSortIndex = m_pnSortIndex;

    nFirst = 0;
    nLast = nGetNumReflns()-1;
    nPass = 0;

    while (1) {
        fFirstValue = (*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].fGetField(nField);
        fLastValue = (*this)[(pnSortIndex)?(pnSortIndex[nLast]):nLast].fGetField(nField);

        if ((nPass==0) && (fLastValue<fValue - fTol)) {
            // If we are finding the first greater, then this reflection gets put at the end of the list.
            if (bFindFirstGreater)
                return m_nNumReflns;
            else
                return -1;
        }

        if ((nPass==0) && (fFirstValue>fValue + fTol)) {
            // If we are finding the first greater, then this reflection gets put at the beginning of the list.
            if (bFindFirstGreater)
                return 0;
            else
                return -1;
        }
        if (nFirst==nLast)
            break;
        
        if ((float)fabs((double)(fFirstValue-fValue)) < fTol)
            break;

        nTest = (nFirst+nLast)/2+1;

        fTestValue = (*this)[(pnSortIndex)?(pnSortIndex[nTest]):nTest].fGetField(nField);
        if (fTestValue>fValue+fTol)
            nLast = nTest-1;
        else 
            nFirst = nTest;

        nPass++;
    }
    if (bFindFirstGreater) {
        // We know that *[nFirst] <= nValue
        // We need to look forward to the first element that is > nValue.
        while ((nFirst+1<m_nNumReflns) && ((float)fabs((double)(*this)[(pnSortIndex)?(pnSortIndex[nFirst+1]):(nFirst+1)].fGetField(nField) -fValue)<=fTol))
            nFirst++;
        if (nFirst<m_nNumReflns)
            nFirst++;
        return nFirst;
    } else {
        
        // Now backtrack to find the first value that is equal to the search value.
        if ((nFirst>=nGetNumReflns()) || ((float)fabs((double)(*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].fGetField(nField)-fValue)>fTol)) 
            return -1;
        while ((nFirst>0) && (((*this)[(pnSortIndex)?(pnSortIndex[nFirst-1]):(nFirst-1)].fGetField(nField)-fValue)<=fTol))
            nFirst--;
        return nFirst;
    }

}

int  Creflnlist::nFindFirst(int nField,const Cstring& sValue,int* pnIndex,bool bFindFirstGreater) {
    int* pnSortIndex;
    int  nFirst,nLast,nTest;
    int  nPass;

    if (nGetNumReflns()==0) {
        if (bFindFirstGreater)
            return 0;
        else
            return -1;
    }

    if (pnIndex)
        pnSortIndex = pnIndex;
    else
        pnSortIndex = m_pnSortIndex;

    nFirst = 0;
    nLast = nGetNumReflns()-1;
    nPass = 0;

    while (1) {
        const Cstring& sFirstValue = (*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].sGetField(nField);
        const Cstring& sLastValue = (*this)[(pnSortIndex)?(pnSortIndex[nLast]):nLast].sGetField(nField);

        if ((nPass==0) && (sLastValue<sValue)) {
            // If we are finding the first greater, then this reflection gets put at the end of the list.
            if (bFindFirstGreater)
                return m_nNumReflns;
            else
                return -1;
        }

        if ((nPass==0) && (sFirstValue>sValue)) {
            // If we are finding the first greater, then this reflection gets put at the beginning of the list.
            if (bFindFirstGreater)
                return 0;
            else
                return -1;
        }
        if (nFirst==nLast)
            break;
        
        if (sFirstValue==sValue)
            break;

        nTest = (nFirst+nLast)/2+1;

        const Cstring& sTestValue = (*this)[(pnSortIndex)?(pnSortIndex[nTest]):nTest].sGetField(nField);
        if (sTestValue>sValue)
            nLast = nTest-1;
        else 
            nFirst = nTest;

        nPass++;
    }
    if (bFindFirstGreater) {
        // We know that *[nFirst] <= nValue.
        // We need to look forward to the first element that is > nValue.
        while ((nFirst+1<m_nNumReflns) && ((*this)[(pnSortIndex)?(pnSortIndex[nFirst+1]):(nFirst+1)].sGetField(nField)<=sValue))
            nFirst++;
        if (nFirst<m_nNumReflns)
            nFirst++;
        return nFirst;
    } else {
        // Now backtrack to find the first value that is equal to the search value.
        if ((nFirst>=nGetNumReflns()) || ((*this)[(pnSortIndex)?(pnSortIndex[nFirst]):nFirst].sGetField(nField)!=sValue)) 
            return -1;
        while ((nFirst>0) && ((*this)[(pnSortIndex)?(pnSortIndex[nFirst-1]):(nFirst-1)].sGetField(nField)==sValue))
            nFirst--;
        return nFirst;
    }
}

int
Creflnlist::nReduce(Cspacegroup &oSpace, const int nAnomFlag,const int nCentricCheck)
{
  // Reduce all reflection indices to asymmetric unit
  // Keep original indices
  // Add packed indices field for sorting and new reduced indices fields
  // Use new indices packing routine
  // Add F+, F- field and centrosymmetric reflection field
  // If nAnomFlag != 0, then use 0,1 in 0th bit of nPackedHKL to flag - or +.

  int i, j;
  int nStat = 0;

  int nPackedHKL;
  int nFplusminus;
  int nCentPhase;

  int a3nReducedHKL[3];

  Crefln *poRefln;

  Cstring asNeededFields[] = {  ms_snPackedHKL,
                                ms_snReducedH,
                                ms_snReducedK,
                                ms_snReducedL,
                                ms_snFplusminusFlag,
                                ms_snCentPhase
                             };

  // Determine which fields are missing ...

  Cstring a6sNewFields[6];
  j = 0;
  for (i = 0; i < 6; i++)
    {
      if (0 > nGetFieldIndex(asNeededFields[i]))
        {
          // Do not have it, so place in request list

          a6sNewFields[j] = asNeededFields[i];
          j++;
        }
    }

  // ... and add them

  if (0 < j) nStat = nExpandRefln(j, a6sNewFields);
  if (nStat != 0) return (-1);

  // Double check

  if (   (0 > m_nFI_nPackedHKL)
      || (0 > m_nFI_nReducedH)
      || (0 > m_nFI_nReducedK)
      || (0 > m_nFI_nReducedL)
      || (0 > m_nFI_nFplusminusFlag)
      || (0 > m_nFI_nCentPhase) )
    {
      cout << "Creflnlist::nReduce - Error expanding reflection field!\n";
      return (-1);
    }

  for (i = 0; i < m_nNumReflns; i++)
    {
      poRefln = *(m_ppoTheReflns+i);
      nCentPhase =  oSpace.nReduceHKL(poRefln, a3nReducedHKL,
                                                      &nFplusminus);
      nPackedHKL = poRefln->nPackHKL(a3nReducedHKL);
      if (0 != nAnomFlag)
        {
          // Keep I+ and I- separated

          // Even if the reflection is centric or absent, 
          // we might be assuming a centric spacegroup during a Laue check.
          // The user can specify a non-centric spacegroup by setting the 
          //    nCentricCheck flag to 0
          // BY DEFAULT nCentricCheck is 1, so nCentPhase has precedence
          //    in the next if.  If 0 == nCentricCheck this overrides
          //                     the centric check

          // If 1==nFplusminus, it's an F+, if -1==nFplusminus it's F-
          // If 0 == nCentPhase, refln is acentric, 
          // If -1 == nCentPhase, refln is systematically absent
          // If 0 < nCentPhase, refln is centric

          if ( (0 == nCentPhase) || (0 == nCentricCheck) )
            {
              // If it's acentric OR we want the I+ and I- of centric
              //     reflections kept separate then ...

              if (1 == nFplusminus)  
                {
                  nPackedHKL = nPackedHKL + 1;
                }
              else if (-1 == nFplusminus)
                {
                  // Negate the hkl

                  a3nReducedHKL[0] = -a3nReducedHKL[0];
                  a3nReducedHKL[1] = -a3nReducedHKL[1];
                  a3nReducedHKL[2] = -a3nReducedHKL[2];
                }
            }
        }
      poRefln->vSetField(m_nFI_nReducedH, a3nReducedHKL[0]);
      poRefln->vSetField(m_nFI_nReducedK, a3nReducedHKL[1]);
      poRefln->vSetField(m_nFI_nReducedL, a3nReducedHKL[2]);
      poRefln->vSetField(m_nFI_nPackedHKL, nPackedHKL);
      poRefln->vSetField(m_nFI_nFplusminusFlag, nFplusminus);
      poRefln->vSetField(m_nFI_nCentPhase, nCentPhase);
    }

  // Sort the results on the nPackedHKL field so that symmetry-related
  //  reflections will be "adjacent" in the sorted list

  if (NULL != m_pnSortIndex)
      delete [] m_pnSortIndex;
  
  m_pnSortIndex = new int [max(1, nGetNumReflns()+1)];
  m_pnSortIndex[nGetNumReflns()] = 0;
  
  if (m_nFI_sBatch>=0)
    vSort2(eReflnField_int_type, m_nFI_nPackedHKL, eReflnField_Cstring_type, m_nFI_sBatch, m_pnSortIndex);
  else
    vSort(eReflnField_int_type, m_nFI_nPackedHKL, m_pnSortIndex);
  
  return (0);
}


void
Creflnlist::vUpdateFieldIndices(void)
{
  // Look at the field names and if any are equal to the member static Cstring
  // field names, then update the member m_nFI_* field index for the field name;
  // Do this in a brute force way.  If the field is not present in the
  // reflection list, then the field index is set to -1 by nGetFieldIndex(...)
  // The nGetFieldIndex call below uses the first letter of the passed
  // field name to determine the type (int, float, Cstring).

  m_nFI_nDetNum         = nGetFieldIndex(ms_snDetNum);
  m_nFI_nNonunfFlag     = nGetFieldIndex(ms_snNonunfFlag);
  m_nFI_nH              = nGetFieldIndex(ms_snH);
  m_nFI_nK              = nGetFieldIndex(ms_snK);
  m_nFI_nL              = nGetFieldIndex(ms_snL);
  m_nFI_fIntensity      = nGetFieldIndex(ms_sfIntensity);
  m_nFI_fSigmaI         = nGetFieldIndex(ms_sfSigmaI);
  m_nFI_fProfitIntensity= nGetFieldIndex(ms_sfProfitIntensity);
  m_nFI_fProfitSigmaI   = nGetFieldIndex(ms_sfProfitSigmaI);
  m_nFI_nPartialFlag    = nGetFieldIndex(ms_snPartialFlag);
  m_nFI_fIntensityPlus  = nGetFieldIndex(ms_sfIntensityPlus);
  m_nFI_fSigmaIPlus     = nGetFieldIndex(ms_sfSigmaIPlus);
  m_nFI_fIntensityMinus = nGetFieldIndex(ms_sfIntensityMinus);
  m_nFI_fSigmaIMinus    = nGetFieldIndex(ms_sfSigmaIMinus);
  m_nFI_fObsPx0         = nGetFieldIndex(ms_sfObsPx0);
  m_nFI_fObsPx1         = nGetFieldIndex(ms_sfObsPx1);
  m_nFI_fObsPxPeak0     = nGetFieldIndex(ms_sfObsPxPeak0);
  m_nFI_fObsPxPeak1     = nGetFieldIndex(ms_sfObsPxPeak1);
  
  m_nFI_fObsSharpness   = nGetFieldIndex(ms_sfObsSharpness);
    
  m_nFI_fObsPeakAreaCount   = nGetFieldIndex(ms_sfObsPeakAreaCount);
  m_nFI_fObsPeakBorderCount = nGetFieldIndex(ms_sfObsPeakBorderCount);
  
  m_nFI_fObsRotMid      = nGetFieldIndex(ms_sfObsRotMid);
  m_nFI_fObsRotEnd      = nGetFieldIndex(ms_sfObsRotEnd);
  m_nFI_fObsRotWidth    = nGetFieldIndex(ms_sfObsRotWidth);
  m_nFI_fObsRotSigma    = nGetFieldIndex(ms_sfObsRotSigma);
  m_nFI_fObsImgMid      = nGetFieldIndex(ms_sfObsImgMid);
  m_nFI_fObsXmm         = nGetFieldIndex(ms_sfObsXmm);
  m_nFI_fObsYmm         = nGetFieldIndex(ms_sfObsYmm);
  m_nFI_fObsZmm         = nGetFieldIndex(ms_sfObsZmm);
  m_nFI_fCalcPx0        = nGetFieldIndex(ms_sfCalcPx0);
  m_nFI_fCalcPx1        = nGetFieldIndex(ms_sfCalcPx1);
  m_nFI_fCalcXmm        = nGetFieldIndex(ms_sfCalcXmm);
  m_nFI_fCalcYmm        = nGetFieldIndex(ms_sfCalcYmm);
  m_nFI_fCalcZmm        = nGetFieldIndex(ms_sfCalcZmm);
  m_nFI_fCalcRotStart   = nGetFieldIndex(ms_sfCalcRotStart);
  m_nFI_fCalcRotMid     = nGetFieldIndex(ms_sfCalcRotMid);
  m_nFI_fCalcRotEnd     = nGetFieldIndex(ms_sfCalcRotEnd);
  m_nFI_fCalcRotWidth   = nGetFieldIndex(ms_sfCalcRotWidth);
  m_nFI_fCalcMosCoeffA  = nGetFieldIndex(ms_sfCalcMosCoeffA);
  m_nFI_fCalcMosCoeffB  = nGetFieldIndex(ms_sfCalcMosCoeffB);
  m_nFI_fPartiality     = nGetFieldIndex(ms_sfPartiality);
  m_nFI_fResolution     = nGetFieldIndex(ms_sfResolution);
  m_nFI_fDetResolution  = nGetFieldIndex(ms_sfDetResolution);
  m_nFI_fRecipCoord0    = nGetFieldIndex(ms_sfRecipCoord0);
  m_nFI_fRecipCoord1    = nGetFieldIndex(ms_sfRecipCoord1);
  m_nFI_fRecipCoord2    = nGetFieldIndex(ms_sfRecipCoord2);
  m_nFI_fRecipCoordD0   = nGetFieldIndex(ms_sfRecipCoordD0);
  m_nFI_fRecipCoordD1   = nGetFieldIndex(ms_sfRecipCoordD1);
  m_nFI_fRecipCoordD2   = nGetFieldIndex(ms_sfRecipCoordD2);
  m_nFI_nReflnNum1      = nGetFieldIndex(ms_snReflnNum1);
  m_nFI_nReflnNum2      = nGetFieldIndex(ms_snReflnNum2);
  m_nFI_nOrignlReflnNum = nGetFieldIndex(ms_snOrignlReflnNum);
  m_nFI_nRefineBatch    = nGetFieldIndex(ms_snRefineBatch);
  m_nFI_nDiffFreq       = nGetFieldIndex(ms_snDiffFreq);
  m_nFI_fFloatH         = nGetFieldIndex(ms_sfFloatH);
  m_nFI_fFloatK         = nGetFieldIndex(ms_sfFloatK);
  m_nFI_fFloatL         = nGetFieldIndex(ms_sfFloatL);
  m_nFI_fHKLResid       = nGetFieldIndex(ms_sfHKLResid);
  m_nFI_fPolarz         = nGetFieldIndex(ms_sfPolarz);
  m_nFI_fLorentz        = nGetFieldIndex(ms_sfLorentz);
  m_nFI_fOblique        = nGetFieldIndex(ms_sfOblique);
  m_nFI_fDeltaPx0       = nGetFieldIndex(ms_sfDeltaPx0);
  m_nFI_fDeltaPx1       = nGetFieldIndex(ms_sfDeltaPx1);
  m_nFI_fDeltaRot       = nGetFieldIndex(ms_sfDeltaRot);
  m_nFI_fBackground     = nGetFieldIndex(ms_sfBackground);
  m_nFI_fBackgroundSigma= nGetFieldIndex(ms_sfBackgroundSigma);
  m_nFI_nTwinID         = nGetFieldIndex(ms_snTwinID);
  m_nFI_fGonio1         = nGetFieldIndex(ms_sfGonio1);
  m_nFI_fGonio2         = nGetFieldIndex(ms_sfGonio2);
  m_nFI_fGonio3         = nGetFieldIndex(ms_sfGonio3);
  m_nFI_fGonio4         = nGetFieldIndex(ms_sfGonio4);
  m_nFI_fGonio5         = nGetFieldIndex(ms_sfGonio5);
  m_nFI_nGonioRotAxis   = nGetFieldIndex(ms_snGonioRotAxis);
  m_nFI_sRejectString   = nGetFieldIndex(ms_ssRejectString);

  m_nFI_fSvec[0]        = nGetFieldIndex(ms_sfSvec0);
  m_nFI_fSvec[1]        = nGetFieldIndex(ms_sfSvec1);
  m_nFI_fSvec[2]        = nGetFieldIndex(ms_sfSvec2);
  m_nFI_fS0vec[0]       = nGetFieldIndex(ms_sfS0vec0);
  m_nFI_fS0vec[1]       = nGetFieldIndex(ms_sfS0vec1);
  m_nFI_fS0vec[2]       = nGetFieldIndex(ms_sfS0vec2);

  m_nFI_fEllipsoidA00   = nGetFieldIndex(ms_sfEllipsoidA00);
  m_nFI_fEllipsoidA01   = nGetFieldIndex(ms_sfEllipsoidA01);
  m_nFI_fEllipsoidA11   = nGetFieldIndex(ms_sfEllipsoidA11);
  m_nFI_fEllipsoidb0    = nGetFieldIndex(ms_sfEllipsoidb0);
  m_nFI_fEllipsoidb1    = nGetFieldIndex(ms_sfEllipsoidb1);
  m_nFI_fEllipsoidc     = nGetFieldIndex(ms_sfEllipsoidc);
  m_nFI_fEllipseMajorAxis       = nGetFieldIndex(ms_sfEllipseMajorAxis);
  m_nFI_fEllipseMinorAxis       = nGetFieldIndex(ms_sfEllipseMinorAxis);
  m_nFI_fEllipseAxisMajorOffset = nGetFieldIndex(ms_sfEllipseAxisMajorOffset);
  
  m_nFI_nSourceRefNum   = nGetFieldIndex(ms_snSourceRefNum);
  m_nFI_nCompeteRefNum  = nGetFieldIndex(ms_snCompeteRefNum);
  m_nFI_nPackedHKL      = nGetFieldIndex(ms_snPackedHKL);
  m_nFI_nReducedH       = nGetFieldIndex(ms_snReducedH);
  m_nFI_nReducedK       = nGetFieldIndex(ms_snReducedK);
  m_nFI_nReducedL       = nGetFieldIndex(ms_snReducedL);
  m_nFI_nFplusminusFlag = nGetFieldIndex(ms_snFplusminusFlag);
  m_nFI_nCentPhase      = nGetFieldIndex(ms_snCentPhase);

  m_nFI_sBatch          = nGetFieldIndex(ms_ssBatch);
  m_nFI_sTag            = nGetFieldIndex(ms_ssTag);
}

int
Creflnlist::nWriteHeader(FILE* pFOut, const char *pcSelectIn,const bool bBinary)
{
  // Write the reflnlist header to the open output stream
  // Return any write errors in *poOut

  int i;
  int nx,ny,nz;

  if (0 == m_nNoWriteHeader)
  {
      // Yes, we do want the header
      
      // Write out only selected fields, so first figure out how
      // many are selected.
      
      char *pcTemp     = m_pcSelect; // Normally m_pcSelect will be all 0's.
      const char *pcIn = pcSelectIn;
      int  nI = 0;
      int  nF = 0;
      int  nS = 0;

      if (pcIn==NULL)
          pcIn = m_pcSelect;
      
      for (i = 0; i < m_nIntReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              nI++;
      }
      for (i = 0; i < m_nFloatReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              nF++;
      }
      for (i = 0; i < m_nCstringReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              nS++;
      }
      
      fprintf(pFOut,"%d %d %d",nI,nF,nS);
      
      // Line 2-?: Write out the text header.  This should not be output
      //           if we are writing in binary mode, since this information is stored
      //           in the binary header.
      
      
      if ((!bBinary) && (m_sEmbeddedHeader.length())) {
          // Count the number of '\n' (or '\r') characters.
          for (nx=0,ny=0,nz=0;nx<m_sEmbeddedHeader.length();nx++) {
              if (m_sEmbeddedHeader.GetAt(nx)=='\n')
                  ny++;
              if (m_sEmbeddedHeader.GetAt(nx)=='\r')
                  nz++;
          }
          nz = max(ny,nz);  // Use maximum of either of these character types.
          fprintf(pFOut," %d\n",nz);
          // Write out the embedded header.
          fprintf(pFOut,"%s",m_sEmbeddedHeader.string());
      } else {
          // No embedded header.
          fprintf(pFOut,"\n");
      }
      
      // Line (?+1)-(?+2): Write out names of fields
      
      pcTemp = m_pcSelect;
      pcIn   = pcSelectIn;

      if (pcIn==NULL)
          pcIn = m_pcSelect;

      for (i = 0; i < m_nIntReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              fprintf(pFOut,"%s\n",(m_psIntNames+i)->string());
      }
      for (i = 0; i < m_nFloatReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              fprintf(pFOut,"%s\n",(m_psFloatNames+i)->string());
      }
      for (i = 0; i < m_nCstringReflnFields; i++) {
          if (*pcTemp++ == *pcIn++)
              fprintf(pFOut,"%s\n",(m_psCstringNames+i)->string());
      }
  }

  // Set output stream to always show decimal point for floating point numbers
  // poOut->setf(ios::showpoint);

  return (ferror(pFOut));
}



int Creflnlist::nReadHeader(const Cstring& sFilename,
                            char** ppcBuffer,
                            int nMaxCharsInBuffer,
                            char** ppcExtraData,
                            tagBinary_Ref_Header& oHeader,
                            int& nBytesOrLinesRead) 
{
    static  Cstring     oLabel("");
    int                 nSkipLines = 0;
    int                 nx = 0, ny = 0;
    char*               pcTok = NULL;

  // Read the reflnlist header from the open input stream
  // Return any read error in *poIn

  // Line 1: Read number of int, float and Cstring fields that are in the file.
  //         Danger! as the following are used to show alloc sizes!

    int                 i = 0;                    // Loop counter


  Cstring sTemp;
  pcTok = strtok(*ppcBuffer,"\n\r\t ");
  if ((!pcTok) || (1!=sscanf(pcTok,"%d",&m_nIntReflnFields)))
      return 1;
  pcTok = strtok(NULL,"\n\r\t ");
  if ((!pcTok) || (1!=sscanf(pcTok,"%d",&m_nFloatReflnFields)))
      return 1;
  pcTok = strtok(NULL,"\n\r\t ");
  if ((!pcTok) || (1!=sscanf(pcTok,"%d",&m_nCstringReflnFields)))
      return 1;
  pcTok = strtok(NULL,"\n\r\t ");
  if (!pcTok)
      return 1;
  if (1==sscanf(pcTok,"%d",&nSkipLines)) {
      sTemp="";
      // Read in the embedded header.
      // The first read simply gets the first '\n'
      for (nx=0;nx<nSkipLines;nx++) {
          pcTok = strtok(NULL,"\n\r");
          if (!pcTok)
              return 1;
          sTemp += pcTok;
          sTemp += "\n";
      }
      m_sEmbeddedHeader = sTemp;
      pcTok = NULL;
  } else {
      nSkipLines = 0;  
  }
  
  m_nTotalFieldsPlus1 = m_nIntReflnFields + m_nFloatReflnFields
      + m_nCstringReflnFields + 1;

  nBytesOrLinesRead = nSkipLines + m_nTotalFieldsPlus1;
  
  // Line 2-..., read in names of fields and allocate space for temp vars
  
  if (NULL != m_psIntNames)     delete [] m_psIntNames;
  if (NULL != m_psFloatNames)   delete [] m_psFloatNames;
  if (NULL != m_psCstringNames) delete [] m_psCstringNames;
  if (NULL != m_pcSelect)       delete [] m_pcSelect;
  
  m_psIntNames     = new Cstring [m_nIntReflnFields];
  m_psFloatNames   = new Cstring [m_nFloatReflnFields];
  m_psCstringNames = new Cstring [m_nCstringReflnFields];
  m_pcSelect       = new char    [m_nTotalFieldsPlus1];
  
  for (i = 0; i < m_nTotalFieldsPlus1; i++)
      m_pcSelect[i] = char(0);
  
  // Now read in names of fields
  
  for (i = 0; i < m_nIntReflnFields;    i++) {
      if (!pcTok)
          pcTok = strtok(NULL,"\n\r");
      if ((!pcTok) || (strchr(pcTok,' ')))
          return 1;
      *(m_psIntNames+i) = pcTok;
      pcTok = NULL;
  }
  for (i = 0; i < m_nFloatReflnFields;  i++) {
      pcTok = strtok(NULL,"\n\r");
      if ((!pcTok) || (strchr(pcTok,' ')))
          return 1;
      *(m_psFloatNames+i) = pcTok;
  }
  for (i = 0; i < m_nCstringReflnFields; i++) {
      pcTok = strtok(NULL,"\n\r");
      if ((!pcTok) || (strchr(pcTok,' ')))
          return 1;
      *(m_psCstringNames+i) = pcTok;
  }

  // Next, try to read in the binary header.
  ny = oHeader.nReadHeader(sFilename,
                           ppcBuffer,
                           ppcExtraData,
                           nMaxCharsInBuffer,
                           nx);
  if (ny<0)
      return -1;
  else if (ny>0)
      return 1;
  else {
      nBytesOrLinesRead = nx;
      return 0;
  }
}

int
Creflnlist::nReindex(const double *pfMatIn) {
    float _pfMatIn[3][3];
    vCopyMat3D(pfMatIn,&_pfMatIn[0][0]);
    return nReindex(&_pfMatIn[0][0]);
}
int
Creflnlist::nReindex(const float *pfMatIn)
{
  // Reindex all the hkls in the reflection list by multiplying the
  // hkl by the input matrix *pfMatIn.

  // Should we check determinant of *pfMatIn?  Use return code to signal
  // bad input matrix.

  int     i;
  for (i = 0; i < m_nNumReflns; i++)
    {
      m_ppoTheReflns[i]->vReindex(pfMatIn);
    }
  return (0);
}

int
Creflnlist::nAddResol(Ccrystal &oCrystal)
{
  // Add the ms_sfResolution (fResolution) field to all reflns in the list

  int   i;
  float a3fHKL[3];
  float a3fX[3];
  float a3x3fBMat[3][3];
  float fTemp;

  (void) nExpandGetField(ms_sfResolution);  // Insert field if missing

  if (!oCrystal.bIsAvailable())
    {
      return (-1);
    }

  // Get the crystal B matrix into local variable

  oCrystal.nCalcBMatrix();
  oCrystal.vGetBMatrix(&a3x3fBMat[0][0]);

  Crefln *poRefln;
  for (i = 0; i < m_nNumReflns; i++)
    {
      poRefln = m_ppoTheReflns[i];

      // Compute reciprocal lattice vector, then take length

      a3fHKL[0] = (float) poRefln->nGetH();
      a3fHKL[1] = (float) poRefln->nGetK();
      a3fHKL[2] = (float) poRefln->nGetL();
      vMulMat3DVec3D(a3x3fBMat, a3fHKL, a3fX);
      fTemp = fLenVec3D(a3fX);
      if (0.0 != fTemp)
        {
          poRefln->vSetField(m_nFI_fResolution, (float)1./fTemp);
        }
      else
        {
          poRefln->vSetField(m_nFI_fResolution, (float)99999.99);
        }
    }
  return (0);
}


int
Creflnlist::nAddDetResol(Cimage_header *poHeader)
{
  // Add the ms_sfDetResolution (detector resolution) field to all reflns in
  // the list

  int        i;
  //int        nStat;
  //float      fTemp;
  Cdetector *poDetector;
  Csource   *poSource;
  float      a3fS0[3];
  float      fWavelength;

  if (!poHeader->bIsAvailable())
    {
      return (-1);
    }

  // Read detector info from header

  // Get the prefix for the 1st detector, if there isn't one, then
  //  the null string will be returned, which is OK

  Cstring sTemp = "";

  (void) poHeader->nGetValue(Cdetector::ms_sDetectorNames, 1, &sTemp);

  poDetector = new Cdetector(*poHeader, sTemp);

  if (!poDetector->bIsAvailable())
    {
      delete poDetector;
      return (-1);
    }

  // Read source info from header

  poSource = new Csource(*poHeader);
  if (!poSource->bIsAvailable())
    {
      delete poSource;
      delete poDetector;
      return (-1);
    }

  poSource->vCalcGetS0(a3fS0);
  fWavelength = poSource->m_poWavelength->fGetWavelength(0);

  // Get min and max resolution on the detector.  This forces a
  // call to the detector vCalcGetDDDN routine which will be needed later.

  float fResMin, fResMax, fResMaxEdge;
  poDetector->nGetResolution(a3fS0, &fResMin, &fResMax, &fResMaxEdge);

  // Insert field for the answer, if missing

  m_nFI_fDetResolution = nExpandGetField(ms_sfDetResolution);

  Crefln *poRefln;
  int   nFI_fPx0, nFI_fPx1;
  float f1, f2;

  // Use observed pixel position if available, otherwise use calc position

  nFI_fPx0 = m_nFI_fCalcPx0;
  nFI_fPx1 = m_nFI_fCalcPx1;

  if (0 < m_nFI_fObsPx0)
    nFI_fPx0 = m_nFI_fObsPx0;
  if (0 < m_nFI_fObsPx1)
    nFI_fPx1 = m_nFI_fObsPx1;

  if (  (0 > nFI_fPx0) || (0 > nFI_fPx1) )
    {
      // Pixel positions not available

      delete poDetector;
      delete poSource;
      return (-1);
    }

  float fReso;
  for (i = 0; i < m_nNumReflns; i++)
    {
      poRefln = m_ppoTheReflns[i];

      f1 = poRefln->fGetField(nFI_fPx0);
      f2 = poRefln->fGetField(nFI_fPx1);

      fReso = poDetector->fCalcGetResolution(f1, f2, a3fS0)
                  * fWavelength;
      poRefln->vSetField(m_nFI_fDetResolution, fReso);
    }
  delete poSource;
  delete poDetector;

  return (0);
}


int
Creflnlist::nDeleteFree(void)
{
  // Delete all the reflections in the free list

  int     i;
  for (i = 0; i < m_nNumFree; i++)
    {
      if (NULL != m_ppoFreeReflns[i])
        {
          delete m_ppoFreeReflns[i];
          m_ppoFreeReflns[i] = NULL;
        }
    }

  if (0 < m_nNumAllocFree)
    {
      delete [] m_ppoFreeReflns;
      m_ppoFreeReflns = NULL;
      m_nNumAllocFree = 0;
    }
  m_nNumFree = 0;

  return (0);
}


int Creflnlist::nOverlapCheck(const int nFI_nHead, const int nFI_nNext, double fRotWidth, double fScanRotStart)
{
    const int   nMaxFields = 3;                 // Maximum number of fields to check overlaps on.
    
    double      afMin[nMaxFields];              // Minimum value of each overlap field over entire data set.
    double      afMax[nMaxFields];              // Maximum value of each overlap field over entire data set.
    
    double      afMaxWidth[nMaxFields];         // Maximum width.
    
    int         anBinSize[nMaxFields];          // Number of bins to generate for each overlap field.
    int         anCumulBinSize[nMaxFields];     // Used for easy indexing of array. [0] dimension changes the fastest.
    
    int         nNumReflns;                     // Total number of reflections in the data set.
    int         nNumOverlaps;                   // Total number of reflections having overlaps.
    
    int         nTotalBins;                     // Total Bins.
    
    int         nPow3Fields;                    // Used when checking bins.  
    
    int         nFI_nIndex;                     // Index bin.
    int*        pnIndex;                        // Sort on Index bin.

    int         nField;
    int         nRef,nx,ny;
    double      f0;

    const int nNumFields = 3;
    int         pnFields[3] = { m_nFI_fObsPx0, m_nFI_fObsPx1, m_nFI_fObsRotMid };
    
    nNumOverlaps = 0;

    nNumReflns   = nGetNumReflns();
    nFI_nIndex = nExpandGetField("nIndex");
    
    // Check the legitimacy of the input fields before anything else.
    // Compute the nPow3Fields value.

    if( pnFields[0] < 0 )
        pnFields[0] = m_nFI_fCalcPx0;
    
    if( pnFields[1] < 0 )
        pnFields[1] = m_nFI_fCalcPx1;
    
    if( pnFields[2] < 0 )
        pnFields[2] = m_nFI_fCalcRotMid;
    
    nPow3Fields = 1;
    for(nField = 0; nField < nNumFields; nField++)
    {
        if( !(pnFields[nField] < m_nFloatReflnFields && 0 <= pnFields[nField]) )
        {
            cout << "ERROR, overlap field " << pnFields[nField] << " is not valid!" << endl;
            return -1;
        }
        
        nPow3Fields *= 3;
    }

    if( 0 == nNumReflns || !(*this)[0].bEllipsoidAvailable() )
        return -1;
  
    // Discover the maximum and minimum values for each value.
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        for(nx=0; nx < nNumFields; nx++)
        {
            f0 = (*this)[nRef].fGetField(pnFields[nx]);
            
            if( nRef == 0 || f0 < afMin[nx] )
                afMin[nx] = f0;
            
            if( nRef == 0 || f0 > afMax[nx] )
                afMax[nx] = f0;
        }

        if( 0 != (*this)[nRef].nConvertToAxialEllipse(&f0, NULL) )
            return -1;
        
        if( nRef == 0 || f0 > afMaxWidth[0] )
            afMaxWidth[0] = f0;
        
        if( nRef == 0 || f0 > afMaxWidth[1] )
            afMaxWidth[1] = f0;
    }
    afMaxWidth[2] = fRotWidth;

    // Figure out how many bins we need to compute.  
    nTotalBins = 1;
    for(nx=0; nx < nNumFields; nx++)
    {
        if( afMaxWidth[nx] <= 0.0 )
            return -1;
        
        ny = (int) ceil((afMax[nx] - afMin[nx])/afMaxWidth[nx]);
        
        anBinSize[nx] = max(1, ny);
        
        if( nx > 0 )
            anCumulBinSize[nx] = anBinSize[nx-1]*anCumulBinSize[nx-1];
        else
            anCumulBinSize[nx] = 1;
        
        nTotalBins *= anBinSize[nx];
    }

    // Index each reflection.  
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        nx = 0;
        for(nField = 0; nField < nNumFields; nField++)
        {        
            f0 = (*this)[nRef].fGetField(pnFields[nField]);
            ny = min(anBinSize[nField]-1,(int) ((f0 - afMin[nField]) / afMaxWidth[nField]));
            nx += ny*anCumulBinSize[nField];
        }
        (*this)[nRef].vSetField(nFI_nIndex, nx);
    }

    // Sort the data on the index.
    pnIndex = new int[nNumReflns];
    vSortIntHeap(nFI_nIndex, pnIndex);  

    // Initialize the head and next fields.  Thus, each reflection is
    // in it's own list of length 1.
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        (*this)[nRef].vSetField(nFI_nHead, nRef);
        (*this)[nRef].vSetField(nFI_nNext,-1);
    }

    // Now, pass through the reflections.
    bool    bFirstTimeSetting_nRef_as_Head = false;
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        bFirstTimeSetting_nRef_as_Head = true;

        int nOffsetCode;
        int nOffsetCodeShift;
        int nCenterCode;
        int nCenterCodeShift;
        int nTestCode;
        int nOffset;
        int nCode;
        int nIndex;
        double a7fEllipsoid1[7];
        double a2fCenter1[2];
        double a7fEllipsoid2[7];
        double a2fCenter2[2];
        int nHead;
        int nTail;

        // Obtain the width and center information of this reflection before hand.
        // This way, we will not have to keep on recomputing it.
        (*this)[nRef].nPutGetEllipse(false,&a7fEllipsoid1[0],&a7fEllipsoid1[4],&a7fEllipsoid1[6],&a2fCenter1[0]);
       
        nCenterCode = (*this)[nRef].nGetField(nFI_nIndex);

        // To check for overlaps, we need to look at neighbors.
        for(nOffsetCode=0; nOffsetCode < nPow3Fields; nOffsetCode++)
        {
            nCenterCodeShift = nCenterCode;
            nOffsetCodeShift = nOffsetCode;
            nTestCode = 0;

            // Build the offseted test code.  We will search for reflections that have this test code.
            for(nField=0; nField < nNumFields; nField++)
            {
                nCode = nCenterCodeShift - (nCenterCodeShift / anBinSize[nField]) * anBinSize[nField];
                
                nCenterCodeShift = nCenterCodeShift / anBinSize[nField];
                
                nOffset = nOffsetCodeShift - (nOffsetCodeShift / 3) * 3;
                
                nOffsetCodeShift /= 3;

                // We should have nOffset set to 0,1, or 2.  (Thus we subtract 1 to get -1,0, or 1)
                // We should have nCode equal to the bin of this field and reflection.

                // Perhaps the test code is out of bounds?
                if( nOffset + nCode - 1 < 0 || nOffset + nCode - 1 >= anBinSize[nField] )
                    break;
                
                nTestCode += (nOffset+nCode-1) * anCumulBinSize[nField];                    
            }

            // Did we make it all the way through all of the fields? If not, then halt.
            if( nField != nNumFields )
                continue;            

            // Do a binary search to obtain the first index having the same code.
            nIndex = nFindFirst(nFI_nIndex, nTestCode, pnIndex);
            if( nIndex >= 0 )
            { 
                // Search all reflections having the same index
                for (;(nIndex<nNumReflns) && ((*this)[pnIndex[nIndex]].nGetField(nFI_nIndex)==nTestCode);nIndex++)
                {
                    // We might be looking at the identical reflection.  Ignore if so.
                    if( pnIndex[nIndex] == nRef )
                        continue;

                    // Do not even check if the 'head' value is the same for both reflections.
                    // In such a case we already found the intersection on a previous step. 
                    if( (*this)[nRef].nGetField(nFI_nHead) == (*this)[pnIndex[nIndex]].nGetField(nFI_nHead) ) 
                        continue;
                    
                    if( ABS((*this)[nRef].fGetField(pnFields[2]) - (*this)[pnIndex[nIndex]].fGetField(pnFields[2])) <= fRotWidth * 1.01 )
                    {
                        // RB 02/22/06 Additional test to make sure the reflections do belong to the same oscillation frame
                        if( -999.0 != fScanRotStart && !bTwoReflnsOnSameFrame(nRef, pnIndex[nIndex], fScanRotStart, fRotWidth) )
                            continue;

                        // See if the two reflections intersect.
                        (*this)[pnIndex[nIndex]].nPutGetEllipse(false, &a7fEllipsoid2[0], &a7fEllipsoid2[4], &a7fEllipsoid2[6], &a2fCenter2[0]);
                        
                        if( (f0 = fEvalEllipse(a2fCenter1[0] - a2fCenter2[0],
                                               a2fCenter1[1] - a2fCenter2[1],
                                               &a7fEllipsoid2[0],
                                               &a7fEllipsoid2[4],
                                               &a7fEllipsoid2[6])) <= 1.0      ||
                            (f0 = fEvalEllipse(a2fCenter2[0] - a2fCenter1[0],
                                               a2fCenter2[1] - a2fCenter1[1],
                                               &a7fEllipsoid1[0],
                                               &a7fEllipsoid1[4],
                                               &a7fEllipsoid1[6])) <= 1.0
                           ) 
                        {
                            
                            // We have found an overlap.  
                            // One of the reflections becomes the head, and the other get's
                            // tacked onto the end of the head's list.  
                            
                            /////////////////////////////////////////////////////////////////////////////////////////
                            // Update the overlap count. First count the reflection with index pnIndex[nIndex],
                            // if we haven't counted it before.
                            if( (*this)[pnIndex[nIndex]].nGetField(nFI_nHead) == pnIndex[nIndex] && // that's what it was initialized to
                                -1 == (*this)[pnIndex[nIndex]].nGetField(nFI_nNext) )  // that's what it was initialized to
                            {
                                nNumOverlaps++;
                            }
                                
                            // We must also count the reflection with index nRef as overlapped,
                            // if we haven't counted it before.
                            if( bFirstTimeSetting_nRef_as_Head && 
                                (*this)[nRef].nGetField(nFI_nHead) == nRef &&  // that's what it was initialized to
                                -1 == (*this)[nRef].nGetField(nFI_nNext) )     // that's what it was initialized to
                            {
                                nNumOverlaps++;
                                
                                bFirstTimeSetting_nRef_as_Head = false;
                            }
                            ////////////////////////////////////////////////////////////////////////////////////////

                            nx = (*this)[nRef].nGetField(nFI_nHead);
                            ny = (*this)[pnIndex[nIndex]].nGetField(nFI_nHead);
                            
                            if( nx < ny )
                            {
                                nHead = nx;
                                nTail = ny;
                            } 
                            else
                            {
                                nHead = ny;
                                nTail = nx;
                            }
                            
                            // Find the end of the 'Head' list.
                            for(nx = nHead; -1 != (ny=(*this)[nx].nGetField(nFI_nNext)); nx = ny);
                            
                            // Set the Next pointer at the end of the 'Head' list.
                            (*this)[nx].vSetField(nFI_nNext, nTail);
                            
                            // Change all reflections in the 'Tail' list to have the correct 'Head' pointer.
                            while( nTail != -1 )
                            {
                                (*this)[nTail].vSetField(nFI_nHead,nHead);
                                nTail = (*this)[nTail].nGetField(nFI_nNext);
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] pnIndex;

#ifdef RB_DEBUG
    
    // Double-check the overlap count
    
    int     nNumOverlapsCheck = 0;

    int     jj = 0;

    for(int ii=0; ii < m_nNumReflns; ii++)
    {
        // If this reflection is the "head" reflection, process it and all the "next" reflections if any.
        if( ii != (*this)[ii].nGetField(nFI_nHead) )
            continue;

        if( -1 != (*this)[ii].nGetField(nFI_nNext) ) // reflection ii overlaps with another reflection
        {
            nNumOverlapsCheck++;

            //go through the linked list and count all overlapped reflections
            jj = ii;
            while( -1 != (jj = (*this)[jj].nGetField(nFI_nNext)) )
            {
                nNumOverlapsCheck++;
            }
        }
    }

#endif

    return nNumOverlaps;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB 10/21/05 The differences between this new function and the previous version of nOverlapCheck() are as follows.
// This function uses an input variable dSpotIntegrBoxWidth instead of reflection ellipsoids. This gives it
// flexibility to do the overlap checking for strategy, when the ellipsoids for all reflections may not be determined yet.
// Also this function makes sure that reflections are physically on the same image before they could be considered
// overlapped.
int Creflnlist::nOverlapCheck(double dSpotIntegrBoxWidth, double fRotWidth, double  fScanRotStart, bool bDisplay)
{
    const int   nMaxFields = 3;                 // Maximum number of fields to check overlaps on.
    
    double      afMin[nMaxFields];              // Minimum value of each overlap field over entire data set.
    double      afMax[nMaxFields];              // Maximum value of each overlap field over entire data set.
    
    double      afMaxWidth[nMaxFields];         // Maximum width.
    
    int         anBinSize[nMaxFields];          // Number of bins to generate for each overlap field.
    int         anCumulBinSize[nMaxFields];     // Used for easy indexing of array. [0] dimension changes the fastest.
    
    int         nNumReflns;                     // Total number of reflections in the data set.
    int         nNumOverlaps;                   // Total number of reflections having overlaps.
    
    int         nTotalBins;                     // Total Bins.
    
    int         nPow3Fields;                    // Used when checking bins.  
    
    int*        pnIndex;                        // Sort on Index bin.

    int         nField;
    int         nRef,nx,ny;
    double      f0;

    const int   nNumFields = 3;
    int         pnFields[3] = { m_nFI_fObsPx0, m_nFI_fObsPx1, m_nFI_fObsRotMid };
    
    nNumOverlaps = 0;
    nNumReflns   = nGetNumReflns();
    
    int     nFI_nHead  = nExpandGetField("nMergeHead");
    int     nFI_nNext  = nExpandGetField("nMergeNext");
    int     nFI_nIndex = nExpandGetField("nIndex");
    
    // Check the legitimacy of the input fields before anything else.
    // Compute the nPow3Fields value.
    if( pnFields[0] < 0 )
        pnFields[0] = m_nFI_fCalcPx0;
    
    if( pnFields[1] < 0 )
        pnFields[1] = m_nFI_fCalcPx1;
    
    if( pnFields[2] < 0 )
        pnFields[2] = m_nFI_fCalcRotMid;
    
    nPow3Fields = 1;
    for(nField = 0; nField < nNumFields; nField++)
    {
        if( !(pnFields[nField] < m_nFloatReflnFields && 0 <= pnFields[nField]) )
        {
            cout << "ERROR, overlap field " << pnFields[nField] << " is not valid!" << endl;
            return -1;
        }
        
        nPow3Fields *= 3;
    }

    if( 0 == nNumReflns )
    {
        printf("/nERROR in reflection overlap checking: no reflections found in the list.\n");
        return -1;
    }
    // Discover the maximum and minimum values for each value.
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        for(nx=0; nx < nNumFields; nx++)
        {
            f0 = (*this)[nRef].fGetField(pnFields[nx]);
            
            if( nRef == 0 || f0 < afMin[nx] )
                afMin[nx] = f0;
            
            if( nRef == 0 || f0 > afMax[nx] )
                afMax[nx] = f0;
        }
    }
    
    /////////////////////////////////////////
    // Assign the widths
    afMaxWidth[0] = dSpotIntegrBoxWidth;
    afMaxWidth[1] = dSpotIntegrBoxWidth;
    afMaxWidth[2] = fRotWidth;
    //////////////////////////////////////////

    // Figure out how many bins we need to compute.  
    nTotalBins = 1;
    for(nx=0; nx < nNumFields; nx++)
    {
        if( afMaxWidth[nx] <= 0.0 )
        {
            printf("ERROR: cannot test reflection overlap, a test difference is zero or negative.\n");
            return -1;
        }
        ny = (int) ceil((afMax[nx] - afMin[nx]) / afMaxWidth[nx]);
        
        anBinSize[nx] = max(1, ny);
        
        if( nx > 0 )
            anCumulBinSize[nx] = anBinSize[nx-1]*anCumulBinSize[nx-1];
        else
            anCumulBinSize[nx] = 1;
        
        nTotalBins *= anBinSize[nx];
    }

    // Index each reflection.  
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        nx = 0;
        for(nField = 0; nField < nNumFields; nField++)
        {        
            f0 = (*this)[nRef].fGetField(pnFields[nField]);
            ny = min(anBinSize[nField]-1,(int) ((f0 - afMin[nField]) / afMaxWidth[nField]));
            nx += ny*anCumulBinSize[nField];
        }
        (*this)[nRef].vSetField(nFI_nIndex, nx);
    }

    // Sort the data on the index.
    pnIndex = new int[nNumReflns];
    vSortIntHeap(nFI_nIndex, pnIndex);  

    // Initialize the head and next fields.  Thus, each reflection is
    // in it's own list of length 1.
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        (*this)[nRef].vSetField(nFI_nHead, nRef);
        (*this)[nRef].vSetField(nFI_nNext,-1);
    }

    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    // Now, pass through the reflections.
    int     nOffsetCode = -1;
    int     nOffsetCodeShift = -1;
    int     nCenterCode = -1;
    int     nCenterCodeShift = -1;
    int     nTestCode = -1;
    int     nOffset = -1;
    int     nCode = -1;
    int     nIndex = -1;
    int     nHead = -1;
    int     nTail = -1;
    bool    bOverlap = false;

    int     nPrintDotStep = max(1, nNumReflns / 49);  // 49, because we want 50 points printed

    printf("\n");
    bool    bFirstTimeSetting_nRef_as_Head = false;
    
    for(nRef=0; nRef < nNumReflns; nRef++)
    {
        if( fmod((double)nRef, (double)nPrintDotStep) == 0.0 )
            printf(".");

        bFirstTimeSetting_nRef_as_Head = true;    
        
        nCenterCode = (*this)[nRef].nGetField(nFI_nIndex);

        // To check for overlaps, we need to look at neighbors.
        for(nOffsetCode=0; nOffsetCode < nPow3Fields; nOffsetCode++)
        {
            nCenterCodeShift = nCenterCode;
            nOffsetCodeShift = nOffsetCode;
            nTestCode = 0;

            // Build the offseted test code.  We will search for reflections that have this test code.
            for(nField=0; nField < nNumFields; nField++)
            {
                nCode = nCenterCodeShift - (nCenterCodeShift / anBinSize[nField]) * anBinSize[nField];
                
                nCenterCodeShift = nCenterCodeShift / anBinSize[nField];
                
                nOffset = nOffsetCodeShift - (nOffsetCodeShift / 3) * 3;
                
                nOffsetCodeShift /= 3;

                // We should have nOffset set to 0,1, or 2.  (Thus we subtract 1 to get -1,0, or 1)
                // We should have nCode equal to the bin of this field and reflection.

                // Perhaps the test code is out of bounds?
                if( nOffset + nCode - 1 < 0 || nOffset + nCode - 1 >= anBinSize[nField] )
                    break;
                
                nTestCode += (nOffset+nCode-1) * anCumulBinSize[nField];                    
            }

            // Did we make it all the way through all of the fields? If not, then halt.
            if( nField != nNumFields )
                continue;            

            // Do a binary search to obtain the first index having the same code.
            nIndex = nFindFirst(nFI_nIndex, nTestCode, pnIndex);
            if( nIndex >= 0 )
            { 
                // Search all reflections having the same index
                for (;(nIndex<nNumReflns) && ((*this)[pnIndex[nIndex]].nGetField(nFI_nIndex)==nTestCode);nIndex++)
                {
                    // We might be looking at the identical reflection.  Ignore if so.
                    if( pnIndex[nIndex] == nRef )
                        continue;

                    // Do not even check if the 'head' value is the same for both reflections.
                    // In such a case we already found the intersection on a previous step. 
                    if( (*this)[nRef].nGetField(nFI_nHead) == (*this)[pnIndex[nIndex]].nGetField(nFI_nHead) ) 
                        continue;
                    
                    bOverlap = true;
                    for(int iField=0; iField < nMaxFields; iField++)
                    {
                        if( ABS((*this)[nRef].fGetField(pnFields[iField]) - 
                            (*this)[pnIndex[nIndex]].fGetField(pnFields[iField])) > afMaxWidth[iField] )
                        {
                            bOverlap = false;
                            break;
                        }    
                    }
                    
                    if( !bOverlap )
                        continue;
                    
                    // RB 02/22/06 Additional test to make sure the reflections do belong to the same oscillation frame
                    if( !bTwoReflnsOnSameFrame(nRef, pnIndex[nIndex], fScanRotStart, fRotWidth) )
                        continue;
                    
                    // We have found an overlap.  
                    // One of the reflections becomes the head, and the other get's
                    // tacked onto the end of the head's list.  

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Update the overlap count. First count the reflection with index pnIndex[nIndex],
                    // if we haven't counted it before.
                    if( (*this)[pnIndex[nIndex]].nGetField(nFI_nHead) == pnIndex[nIndex] && // that's what it was initialized to
                        -1 == (*this)[pnIndex[nIndex]].nGetField(nFI_nNext) )  // that's what it was initialized to
                    {
                        nNumOverlaps++;
                    }
                        
                    // We must also count the reflection with index nRef as overlapped,
                    // if we haven't counted it before.
                    if( bFirstTimeSetting_nRef_as_Head && 
                        (*this)[nRef].nGetField(nFI_nHead) == nRef &&  // that's what it was initialized to
                        -1 == (*this)[nRef].nGetField(nFI_nNext) )     // that's what it was initialized to
                    {
                        nNumOverlaps++;
                        
                        bFirstTimeSetting_nRef_as_Head = false;
                    }
                    ////////////////////////////////////////////////////////////////////////////////////////

                    nx = (*this)[nRef].nGetField(nFI_nHead);
                    ny = (*this)[pnIndex[nIndex]].nGetField(nFI_nHead);
                    if( nx < ny )
                    {
                        nHead = nx;
                        nTail = ny;
                    } 
                    else
                    {
                        nHead = ny;
                        nTail = nx;
                    }

                    // Find the end of the 'Head' list.
                    for(nx = nHead; -1 != (ny=(*this)[nx].nGetField(nFI_nNext)); nx = ny);

                    // Set the Next pointer at the end of the 'Head' list.
                    (*this)[nx].vSetField(nFI_nNext, nTail);

                    // Change all reflections in the 'Tail' list to have the correct 'Head' pointer.
                    while( nTail != -1 )
                    {
                        (*this)[nTail].vSetField(nFI_nHead,nHead);
                        nTail = (*this)[nTail].nGetField(nFI_nNext);
                    }
                }
            }
        }
    }
    printf("\n");

    delete[] pnIndex;
    

    if( nNumOverlaps > 0 && bDisplay )
    {
        // Go back and set the nonunf field so that the display program could display the overlapped reflections.
        int     nNonunfField  = nExpandGetField(ms_snNonunfFlag);
        for(nx = 0; nx < nGetNumReflns(); nx++)
        {
            if( (*this)[nx].nGetField(nFI_nNext) != -1 || (*this)[nx].nGetField(nFI_nHead) != nx  )
            {
                (*this)[nx].vSetField(nNonunfField, enColorMagenta);
            }
        }
    }

    return nNumOverlaps;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Creflnlist::LoadViewOptions(Cstring csInColor, Cstring csBeforeColor,
                        Cstring csAfterColor, Cstring csObsShape, Cstring csCalShape)
{
        m_csInColor = csInColor;
        m_csBeforeColor = csBeforeColor;
        m_csAfterColor = csAfterColor;
        m_csObsShape = csObsShape;
        m_csCalShape = csCalShape;
}

int
Creflnlist::nDeleteSpot(float fPx0, float fPx1)
{
  // Delete the refln closest to pixel position (fPx0, fPx1)
  // Returns the number of the refln deleted, or -1 on error

  Crefln *poCurRef;
  int i;
  int nClosestIndex = -1;
  double dCurX, dCurY;
  double dMinDistSqrd = 625.0;
  double dCurDistSqrd;
  int nXIndex, nYIndex;

  if (-1 != m_nFI_fObsPx0)
    {
      nXIndex = m_nFI_fObsPx0;
      nYIndex = m_nFI_fObsPx1;
    }
  else
    {
      nXIndex = m_nFI_fCalcPx0;
      nYIndex = m_nFI_fCalcPx1;
    }

  for (i = 0; i < m_nNumReflns; ++i)
    {
      poCurRef     = m_ppoTheReflns[i];
      dCurX        = poCurRef->fGetField(nXIndex);
      dCurY        = poCurRef->fGetField(nYIndex);
      dCurDistSqrd =  (fPx0 - dCurX) * (fPx0 - dCurX)
                    + (fPx1-dCurY) * (fPx1 - dCurY);
      if (dCurDistSqrd < dMinDistSqrd)
        {
          nClosestIndex = i;
          dMinDistSqrd  = dCurDistSqrd;
        }
    }

  if (-1 != nClosestIndex)
    {
      // Shouldn't it also depend on dMinDistSqrd? That is
      // dMinDistSqrd <= dValue?
      //  in a sense it does: no deletion occurs if dMinDistSqrd >625

      nDelete(nClosestIndex);
      m_bNeedsSaving = 1;
    }
  return (nClosestIndex);
}

void
Creflnlist::vAddSpot(float fPx0, float fPx1, float fIntensity, float fSigmaI,
                     float fReso, float fRotMid, float fRotWidth)
{
  Crefln *poRefln;
  int nXIndex, nYIndex;

  if (-1 < m_nFI_fObsPx0)
    {
      nXIndex = m_nFI_fObsPx0;
      nYIndex = m_nFI_fObsPx1;
    }
  else if (-1 < m_nFI_fCalcPx0)
    {
      nXIndex = m_nFI_fCalcPx0;
      nYIndex = m_nFI_fCalcPx1;
    }
  else
    {
      return;
    }

  // Make a copy of one of the reflections in the list
  // to get some default settings and then fill in the new
  // values for the new reflection.

  poRefln = new Crefln(this);
  poRefln->vSetIntensity(fIntensity);
  poRefln->vSetSigmaI(fSigmaI);
  poRefln->vSetField(nXIndex, fPx0);
  poRefln->vSetField(nYIndex, fPx1);
  if (-1 < m_nFI_fDetResolution)
    {
      poRefln->vSetField(m_nFI_fDetResolution, fReso);
    }
  if (-1 < m_nFI_fObsRotMid)
    {
      poRefln->vSetField(m_nFI_fObsRotMid, fRotMid);
    }
  else if (-1 < m_nFI_fCalcRotMid)
    {
      poRefln->vSetField(m_nFI_fCalcRotMid, fRotMid);
    }

  if (-1 < m_nFI_fObsRotWidth)
    {
      poRefln->vSetField(m_nFI_fObsRotWidth, fRotWidth);
    }
  else if (-1 < m_nFI_fCalcRotWidth)
    {
      poRefln->vSetField(m_nFI_fCalcRotWidth, fRotWidth);
      if (-1 < m_nFI_fCalcRotStart)
        poRefln->vSetField(m_nFI_fCalcRotStart, (float)(fRotMid - 0.5 * fRotWidth));
      if (-1 < m_nFI_fCalcRotEnd)
        poRefln->vSetField(m_nFI_fCalcRotEnd, (float)(fRotMid + 0.5 * fRotWidth));
    }
        
  // Now insert the new reflection at the beginning of
  // the list.  This call will expand the list and
  // increment the number of reflections and will
  // allocate memory for a new reflection so we need to
  // free the memory that we allocated above for the new
  // reflection.

  nInsert(poRefln, 0);
  delete poRefln;
  m_bNeedsSaving = 1;
}


Crefln* Creflnlist::poGetClosestRefln(float x, float y,int nCoord1,int nCoord2)
{
  Crefln *CurRef;
  Crefln *retRef = NULL;
  int i;
  int XIndex, YIndex;
  int ClosestIndex = -1;
  double CurX, CurY;
  double MinDistSqrd = 9999.9;
  double CurDistSqrd;

  if ((nCoord1!=-1) || (nCoord2!=-1))
    {
      XIndex = nCoord1;
      YIndex = nCoord2;
    } 
  else if(m_nFI_fObsPx0 != -1)
    {
      XIndex = m_nFI_fObsPx0;
      YIndex = m_nFI_fObsPx1;
    }
  else
    {
      XIndex = m_nFI_fCalcPx0;
      YIndex = m_nFI_fCalcPx1;
    }
  
  for(i=0; i<m_nNumReflns; ++i)
    {
      CurRef = m_ppoTheReflns[i];
      CurX = CurRef->fGetField(XIndex);      
      CurDistSqrd = (x-CurX)*(x-CurX);
      if (YIndex!=-1) {
        CurY = CurRef->fGetField(YIndex);
        CurDistSqrd += (y-CurY)*(y-CurY);
      }
      if(CurDistSqrd < MinDistSqrd)
        {
          ClosestIndex = i;
          MinDistSqrd = CurDistSqrd;
        }
    }
  if(ClosestIndex != -1)
    {
      retRef = m_ppoTheReflns[ClosestIndex];
    }
  return retRef;
}

bool Creflnlist::bNeedsSaving()
{
  return m_bNeedsSaving;
}
void Creflnlist::NeedsSaving(bool NewVal)
{
  m_bNeedsSaving = NewVal;
}

/*  The diffs function assumes that we already have the sTag flag filled with an appropriately formated tag.

*/
                    

int Creflnlist::nDiffs(int nNumConstantFields,
                       eReflnFieldType* peConstantFields,
                       int* pnConstantFields,
                       float* pfConstantFieldTolerances)
{
    Cspacegroup oSpace;
    Cstring sTemp,sTemp2;
    float f0,f1;

    // HKL Equivalency loop.
    int nLastHKL,nThisHKL;
    int nStartIndex,nEndIndex,nGroupSize;
    int nRef,nRef2,nRef1,nRefSort;
    int nx;
    int nPackedHKL;
    int* pnIndex;
    int nNumRefs;
    int nTagSize;
    int nConstant;



    if (m_nFI_sTag==-1) 
        return 1;

    sTemp = (*this)[0].sGetField(m_nFI_sTag);
    if (!strchr(sTemp.string(),'_'))
        return 1;
    nTagSize=strlen(sTemp.string())-strlen(strchr(sTemp.string(),'_'));

    // Sort on asymmetric unit.  Use spacegroup 1 with I+!=I-
    oSpace.vSet(1);  
    printf("Sorting HKL's into unique groups...");
    nReduce(oSpace,1,0);
    printf("done.\n");
   
    pnIndex = pnGetSortIndex();
    nNumRefs = nGetNumReflns();

    nPackedHKL = m_nFI_nPackedHKL;
    if (nPackedHKL<0)
        return 1;

    nLastHKL = (*this)[pnIndex[0]].nGetField(nPackedHKL);
    nStartIndex = 0;

    for (nRefSort=0;nRefSort<=nNumRefs;nRefSort++)  {
            nRef=pnIndex[nRefSort];
                
            if (nRefSort<nNumRefs)
                nThisHKL = (*this)[nRef].nGetField(nPackedHKL);
            if ((nRefSort!=nNumRefs) && (nThisHKL==nLastHKL))
                continue;
            nEndIndex=nRefSort-1;
            nGroupSize = nEndIndex - nStartIndex + 1;

            float fIntensity = (rand() % 1000 ) + 500;
            
            // Traverse through each of the reflections in this equivalency group, and determine which file ID's they
            // have in common.

            for (nRef1=nStartIndex;nRef1<=nEndIndex;nRef1++) {
                sTemp = (*this)[pnIndex[nRef1]].sGetField(m_nFI_sTag);

                // Now look at the other reflections in the group.
                for (nRef2=nStartIndex;nRef2<=nEndIndex;nRef2++) {
                    if (nRef1!=nRef2) {
                        // Look at each index that is required for equivalence.
                        for (nConstant=0;nConstant<nNumConstantFields;nConstant++) {
                            if (peConstantFields[nConstant]==eReflnField_int_type) {
                                if ((*this)[pnIndex[nRef1]].nGetField(pnConstantFields[nConstant]) !=
                                    (*this)[pnIndex[nRef2]].nGetField(pnConstantFields[nConstant]))
                                    break;
                            } else if (peConstantFields[nConstant]==eReflnField_Cstring_type) {
                                if (!((*this)[pnIndex[nRef1]].sGetField(pnConstantFields[nConstant]) ==
                                    (*this)[pnIndex[nRef2]].sGetField(pnConstantFields[nConstant])))
                                    break;
                            } else if (peConstantFields[nConstant]==eReflnField_float_type) {
                                f0 = (*this)[pnIndex[nRef1]].fGetField(pnConstantFields[nConstant]);
                                f1 = (*this)[pnIndex[nRef2]].fGetField(pnConstantFields[nConstant]);
                                if (fabs((double)(f0-f1)) > pfConstantFieldTolerances[nConstant]*fabs((double)(f0+f1)))
                                    break;
                            } else 
                                return 1;
                        }
                        if (nConstant==nNumConstantFields) {
                            // We have an equivalence.  Set the appropriate string elements.
                            sTemp2 = (*this)[pnIndex[nRef2]].sGetField(m_nFI_sTag);
                            for (nx=0;nx<nTagSize;nx++) {
                                if (sTemp2.string()[nx]=='1')
                                    sTemp.string()[nx]='1';
                            }
                        }
                    }
                }
                // Set the string again.
                (*this)[pnIndex[nRef1]].vSetField(m_nFI_sTag,sTemp);
            }
       
            nStartIndex=nRefSort;
            if (nRefSort<nNumRefs)
                nLastHKL =  (*this)[nRef].nGetField(nPackedHKL);
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



/* Format of log file:

  Each line has the following fields.
  nLevel [ 'FIELDS' sFileName nNumReflns nIntegers nFloats nStrings [ sField1 sField2 ...]] [ 'PROG' nNumArgs [ sProg1 sProg2 ...] ] [['REM' nNumArgs [ sComment1 sComment2 ...] ] [ Further 'REM' statements] ... ]

  nLevel:       The level.  For a reflection list, Level 1 contains data for the reflection list from which this reflection list was derived.
                Level 0 contains data for the reflection list.  
  nNumReflns:   The number of reflections in the list at that point.
  nIntegers,nFloats,nStrings,sField1..sFieldN:  
                Fields matching those in the text header.
  sProg1...sProgN:  
                Command line arguments when the program was written.
  sComment1..sCommentN:  
                Comments
  
*/

int Creflnlist::nDisplayLog(int nVerbose) {
    int nx,ny;
    int nTab;
    static Cstring sBuf;
    static Cstring sOutBuf;
    char* pcTok;
    bool bTerminationOkay;

#define DISPTOK { pcTok = strtok(NULL,"\t "); if (!pcTok) goto loop_terminate; }

    cout << flush;
    cout << flush;
    printf("\nReflection Log\n"
           "--------------\n");

    sBuf = m_sChildLog;
    bTerminationOkay = TRUE;
    for (pcTok = strtok(sBuf.string(),"\t ");pcTok;) {
        
        bTerminationOkay = FALSE;
        sOutBuf = "";

        nTab = atoi(pcTok)*4;
        // Don't let the tab's get out of control!
        nTab = min(4*6,nTab);
        if ((nTab<0) || (nTab>30))
            break;
        DISPTOK;        
        if (!strcmp(pcTok,"FIELDS")) {
            DISPTOK;
            sOutBuf += "Name:    ";
            sOutBuf += pcTok;
            DISPTOK;
            sOutBuf += " ( ";
            sOutBuf += atoi(pcTok);
            sOutBuf += " Reflns) ";
            sOutBuf += "( ";
            for (nx=0,ny=0;nx<3;nx++) {
                DISPTOK;
                ny += atoi(pcTok);
                sOutBuf += atoi(pcTok);
                sOutBuf += " ";
            }
            sOutBuf += ") ";

            if (nVerbose>=2)                          
                sOutBuf += " \nFields:  ";

            for (nx=0;nx<ny;nx++) {
                DISPTOK;
                if (nVerbose>=2) {
                    sOutBuf += pcTok;
                    sOutBuf += " ";
                }
            }
            DISPTOK;
        } else {
            sOutBuf += "Name: ??? ";
        }
        if (!strcmp(pcTok,"PROG")) {
            if ((nVerbose>=1) && (sOutBuf.string()[sOutBuf.length()-1]!='\n'))
                sOutBuf += " \n";
            sOutBuf += "Program: ";
            DISPTOK;
            ny= atoi(pcTok);
            for (nx=0;nx<ny;nx++) {
                DISPTOK;
                if ((nx==0) || (nVerbose>=1)) {
                    sOutBuf += pcTok;
                    sOutBuf += " ";
                }
            }
            if (nVerbose>=1)
                sOutBuf += " \n";
            DISPTOK;
        } else {
            if ((nVerbose>=1) && (sOutBuf.string()[sOutBuf.length()-1]!='\n'))
                sOutBuf += " \n";
            sOutBuf += "Program: ??? ";
        }

        // There might be more than one comment.
        while (!strcmp(pcTok,"REM")) {
            if (sOutBuf.string()[sOutBuf.length()-1]!='\n')
                sOutBuf += " \n";
            sOutBuf += "Remarks: ";
            DISPTOK;
            ny= atoi(pcTok);
            for (nx=0;nx<ny;nx++) {
                DISPTOK;
                sOutBuf += pcTok;
                sOutBuf += " ";
            }
            DISPTOK;
        }

        cout << sTabFormat(sOutBuf,nTab) << "\n" << flush;

        // There might be other arguments types on the same line,
        // that are not yet supported.
        while (strcmp(pcTok,"\n"))
            DISPTOK;
        bTerminationOkay = TRUE;        
        DISPTOK;
    }
loop_terminate:   
    printf("\nEnd Log\n-------\n");
    if (!bTerminationOkay) {
        printf("Bad log format!\n");
        return 1;
    } else
        return 0;

}

///////////////////////////////////////////////
///////////////////////////////////////////////

// Erase the extra data fields.

int Creflnlist::nEraseExtra() {
    m_sChildLog="";
    m_sCommentLog="";
    m_sEmbeddedHeader="";
    return 0;
}

// Parse the log string comming in from a binary reflection file.
// Get the "extra" information.

int Creflnlist::nParseExtra(char* pcInputLog) {
    m_sChildLog=pcInputLog;
    m_sEmbeddedHeader=(pcInputLog+m_sChildLog.length()+1);
    return 0;
}


// Add extra data from one reflection list to another.

int Creflnlist::nAddExtra(Creflnlist& oList) {
    m_sChildLog+=oList.m_sChildLog;
    m_sCommentLog+=oList.m_sCommentLog;
    if (m_sEmbeddedHeader.length()==0)
        m_sEmbeddedHeader=oList.m_sEmbeddedHeader;
    else if (oList.m_sEmbeddedHeader.length()==0)
        ;   // Do nothing.  The inserted list does not have an embedded header.
    else {
        // We need to make sure that the two reflection lists are compatible.
        // If not, then we should warn the user.
        Ccrystal* poCrystal1;
        Ccrystal* poCrystal2;
        poCrystal1=poGetCrystal();
        poCrystal2=oList.poGetCrystal();
        if ((NULL != poCrystal1) && (NULL != poCrystal2))
        {
            int nx,ny;
            float a6fCell1[6];
            float a6fCell2[6];
            float a6fCellD[6];
            float f0;
            
            // See if the spacegroups are the same.
            nx=poCrystal1->m_poSpacegroup->nGet();
            ny=poCrystal2->m_poSpacegroup->nGet();
            if (nx!=ny) {
                cout << "WARNING:  Spacegroups " << nx << " and  " << ny << " do not match!\n" << flush;
            }
            // See if the cells are the same.
            poCrystal1->vGetCell(a6fCell1);
            poCrystal2->vGetCell(a6fCell2);
            vSubVecNDVecND(6,a6fCell1,a6fCell2,a6fCellD);
            f0 = max(fLenVecND(6,a6fCell1),fLenVecND(6,a6fCell2));
            if ((f0 != 0.0) && (fLenVecND(6,a6fCellD)/f0>0.05)) {
                cout << "WARNING:  Cell parameters differ by more than 5 percent\n" << flush;
            }
            // See if the orrientation angles are the same.
            poCrystal1-> vGetOrientAngles(a6fCell1);
            poCrystal2-> vGetOrientAngles(a6fCell2);
            vSubVec3DVec3D(a6fCell1,a6fCell2,a6fCellD);
            f0 = max(fLenVec3D(a6fCell1),fLenVec3D(a6fCell2));
            if ((f0 != 0.0) && (fLenVec3D(a6fCellD)/f0>0.05)) {
                cout << "WARNING:  Cell orientation angles differ by more than 5 percent\n" << flush;
            }            
        }
        // Now, merge the headers.

        if (!poCrystal1)
            m_sEmbeddedHeader=oList.m_sEmbeddedHeader;

        if (poCrystal1 != NULL)
            delete poCrystal1;
        if (poCrystal2 != NULL)
            delete poCrystal2;
    }
    return 0;
}


///////////////////////////////////////////////
///////////////////////////////////////////////

Ccrystal* Creflnlist::poGetCrystal() {
    Cimage_header oHead;
    Ccrystal* poCrystal;
    int nStat;

    if (m_sEmbeddedHeader.length() && (-1 != m_sEmbeddedHeader.find(D_K_CrystalPrefix)))
        nStat=oHead.nParse(m_sEmbeddedHeader);
    else
        return NULL;
    if (nStat)
        return NULL;
    poCrystal = new Ccrystal(oHead);
    if (poCrystal->bIsAvailable()) {
        return  poCrystal;
    } else {
        delete poCrystal;
        return NULL;
    }
}

Cimage_header* Creflnlist::poGetHeader() {
    Cimage_header* poHeader;
    int nStat;

    poHeader = new Cimage_header();
    if (m_sEmbeddedHeader.length()) {
        nStat=poHeader->nParse(m_sEmbeddedHeader);        
    } else
        nStat = 0;
    if (nStat) {
        delete poHeader;
        return NULL;
    }
    return poHeader;
}

int Creflnlist::nPutHeaderInfo(Cimage_header& oHeader,Cstring sPattern,bool bDeleteOld)
{
    Cimage_header*                  poHeader;
    Cstring                         sTemp;
    std::vector<Cstring>            asKeywords;
    int                             nx = 0;
    
    // Load the existing header.
    poHeader = poGetHeader();
    if (!poHeader)
        return 1;
    // Delete all previous keywords in header.
    if (bDeleteOld)
        poHeader->nDeleteMask(sPattern);

    // Copy keywords over.
    poHeader->nCopyMask(oHeader,sPattern);

    if (poHeader->nFindKeywordMask("*",asKeywords))
    {
        delete poHeader;
        return 1;
    }
    
    m_sEmbeddedHeader = "";
    for (nx = 0; nx < asKeywords.size(); nx++)
    {
        poHeader->nGetValue(asKeywords[nx],&sTemp);
        if (!sTemp.after('\n').length()) {
            m_sEmbeddedHeader += asKeywords[nx];
            m_sEmbeddedHeader += "=";
            m_sEmbeddedHeader += sTemp;
            m_sEmbeddedHeader += ";\n";
        }
    }
    
    delete poHeader;
    return 0;
}

int Creflnlist::nPutCrystal(Ccrystal& oCrystal) {
    Cimage_header oHead;

    if (oCrystal.fCalcVolume()==0.0) {
        // The dark crystal.  
        // Do not place it in embedded header.
        return 0;
    }
    if (oCrystal.nUpdateHeader(&oHead))
        return 1;
    if (nPutHeaderInfo(oHead,D_K_HeaderCrystalInfo,false))
        return 1;

    oCrystal.nUpdateHeader(&oHead);
    
   

    return 0;
}




///////////////////////////////////////////////
///////////////////////////////////////////////



// Function updates the "comment" field in the log file.

int Creflnlist::nLogComment(const Cstring& sLog,bool bAdd) {
    if ((bAdd) && (!m_sCommentLog.length()))
        m_sCommentLog += " \n ";
    else
        m_sCommentLog = "";
    m_sCommentLog += sLog;
    m_sCommentLog += " \n ";
    return 0;
}

// Generate the log string to write into the header.

int Creflnlist::nGenerateLog(Cstring& sLogbuf,const Cstring& sFilename,tagBinary_Ref_Header& oBinaryHeader) {
    int nx;
    int nCount;
    char* pcBuf;
    char* pcTok,*pcTok2;
    bool bAppendNewLog;

    static Cstring sBuf;
    static Cstring sTemp;

    // Should we append new log information, or
    // are we working in "silent" mode.  This allows the user
    // to insert changes without affecting the log file.

    if (m_eWriteLog==eLogoutput_log)
        bAppendNewLog = TRUE;
    else if (m_eWriteLog==eLogoutput_nolog)
        bAppendNewLog = FALSE;
    else if (sGetEnv("DTREK_REFLN_NOLOG").length()!=0) 
        bAppendNewLog = FALSE;
    else
        bAppendNewLog = TRUE;


    if (bAppendNewLog) {
        sLogbuf="0 FIELDS ";
        sLogbuf += sFilename;
        sLogbuf += " ";
        sLogbuf += m_nNumReflns;
        sLogbuf += " ";
        for (nx=0;nx<3;nx++) {
            sLogbuf += oBinaryHeader.a3nFields[nx];
            sLogbuf += " ";
        }
        sLogbuf += oBinaryHeader.ms_sFields;
        
        
        // Paste in the command line arguments. argv[0] will get a little extra attention.
        // We want to remove the path, and remove the .exe if we are on the PC.
        
        if (g_sCommandLine.length()) {
            sLogbuf += "PROG ";
            
            // Count number of arguments in list, and place this as the first argument.
            sBuf = g_sCommandLine;
            for (nx=0,pcTok=strtok(sBuf.string(),"\t \n");pcTok;pcTok=strtok(NULL,"\t \n")) 
                nx++;
            sLogbuf += nx;
            sLogbuf += " ";
            
            // Deal with argv[0]
            sBuf = g_sCommandLine;
            pcBuf = &sBuf.string()[0];
            strtok(pcBuf,"\t ");
            for (pcTok2=strtok(pcBuf,"\\/");pcTok2;pcTok2=strtok(NULL,"\\/"))
                pcTok=pcTok2;
            pcTok = strtok(pcTok,".");         
            sLogbuf += pcTok;
            sLogbuf += " ";
            
            
            // Deal with rest of arguments.
            sBuf = g_sCommandLine;
            pcBuf = &sBuf.string()[0];       
            for (pcTok=strtok(pcBuf,"\t \n"),pcTok=strtok(NULL,"\t \n");pcTok;pcTok=strtok(NULL,"\t \n")) {
                sLogbuf += pcTok;
                sLogbuf +=" ";
            }
        }
        if (m_sCommentLog.length()) {
            
            // Count number of arguments in list, and place this as the first argument.
            sBuf = m_sCommentLog;        
            sTemp = "";
            for (nx=0,pcTok=strtok(sBuf.string(),"\t ");pcTok;pcTok=strtok(NULL,"\t ")) {
                if (!strcmp(pcTok,"\n")) {
                    if ((nx) && (sTemp.length())) {
                        sLogbuf += " REM ";
                        sLogbuf += nx;
                        sLogbuf += " ";
                        sLogbuf += sTemp;
                    }
                    sTemp = "";
                    nx=0;
                } else {
                    sTemp += pcTok;
                    sTemp += " ";
                    nx++;
                }
            }
            if (nx) {
                sLogbuf += " REM ";
                sLogbuf += nx;
                sLogbuf += " ";
                sLogbuf += sTemp;
            }
            
        }
        sLogbuf += " \n ";    // Insert a new line to indicate the end of this log field.
    }

    // Next, we parse the log fields that were inside the reflection list.
    // These need to have their "level" value incremented.
    sBuf = m_sChildLog;
    for (nCount=0,pcTok = strtok(sBuf.string(),"\t ");pcTok;pcTok=strtok(NULL,"\t ")) {
        if (nCount==0) {
            // We expect an integer.
            nx = atoi(pcTok);
            // If we appended new log information, then we must increment each
            // level parameter.  Otherwise, we just copy the same level parameter over.
            // Really, this whole loop could be avoided if bAppendNewLog is FALSE... However,
            // I do it just so the data is processed.
            if (bAppendNewLog)
                sLogbuf += (nx+1);
            else
                sLogbuf += (nx);
            sLogbuf += " ";
            nCount++;
        } else if (strcmp(pcTok,"\n")) {
            sLogbuf += pcTok;
            sLogbuf += " ";
            nCount++;
        } else {
            sLogbuf += " \n ";
            nCount=0;
        }         
    }



    return 0;
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// The \0 is for octal
char* tagBinary_Ref_Header::ms_pcId = "\032\012\014\012\004D*TREK_BIN";

CReffloat tagBinary_Ref_Header::ms_fCodeFloat = 1234.2f;
int tagBinary_Ref_Header::ms_nCodeInt = 1234;

bool tagBinary_Ref_Header::ms_bVax2IEEE;
bool tagBinary_Ref_Header::ms_bSwapFloat;
bool tagBinary_Ref_Header::ms_bSwapLong;
Cstring tagBinary_Ref_Header::ms_sFields;



int tagBinary_Ref_Header::nReadHeader(const Cstring& sFilename,
                                      char** ppcBuffer,
                                      char** ppcExtra,
                                      int nMaxSearch,
                                      int& nCharsRead)
{
    int     nBufferSizeIn = nMaxSearch; // save the buffer size
    char*   pBufferIn = *ppcBuffer;

    int     nVerbose = 2;
    static  Cstring sTemp;

    int     nCount;

    char    cChar;
    int     nx;
    int     nPos;

    if (ppcExtra)
        *ppcExtra=NULL;

    for (nCount=0,nPos=0;(nMaxSearch>0) && (ms_pcId[nCount]!=0);nMaxSearch--,nPos++)
    {
        cChar = pBufferIn[nPos];
        if (cChar==ms_pcId[nCount])
            nCount++;
        else
            nCount=0;
    }
    
    if (ms_pcId[nCount])
    {
        return -1;
    }
    memcpy(this,&pBufferIn[nPos],sizeof(tagBinary_Ref_Header));
    nCharsRead = nPos + sizeof(tagBinary_Ref_Header);

    int nStat=0;

    // Set the flags for reading.  These will indicate which conversions need to be done.

    if (((int) cSizeFloat)==8) {
        double fCode;
        ms_bVax2IEEE=FALSE;
        ms_bSwapFloat=FALSE;
        fCode = *((double*) & _unionbrace._structbrace.a16cSampleFloat[0]);
        if (fCode!=(double) ms_fCodeFloat ) {
            cout << "ERROR:  Did not find floating point conversion for 8 byte floats.\n";
            nStat++;
        }
    } else {
        float fCode;
        float fNewCode;
       
        fCode = fCode = *((float*) & _unionbrace._structbrace.a16cSampleFloat[0]);
        if (fCode!=(float) ms_fCodeFloat) {
            fNewCode = fSwapFloat(fCode);
            if (fNewCode!= (float) ms_fCodeFloat) {
                fNewCode = fVAXtoIEEE(&fCode);
                if (fNewCode != (float) ms_fCodeFloat) {
                    cout << "ERROR:  Did not find floating point conversion for 4 byte floats.\n";
                    nStat++;
                } else {
                    ms_bVax2IEEE=TRUE;
                    ms_bSwapFloat=FALSE;
                }
            } else {
                ms_bVax2IEEE=FALSE;
                ms_bSwapFloat=TRUE;
            }
        } else {
            ms_bVax2IEEE = FALSE;
            ms_bSwapFloat = FALSE;
        }
    }

   
    int nCode;

    nCode = *((int*) & _unionbrace._structbrace.a16cSampleInt[0]);
    if (nCode != ms_nCodeInt) {
        if (nSwapLong(nCode) != ms_nCodeInt) {
            cout << "ERROR:  Did not find conversion for 4 byte integers.\n";
            nStat++;
        } else {
            ms_bSwapLong = TRUE;
        }
    } else {
        ms_bSwapLong = FALSE;
    }

    if (nVerbose > 3 ) {
        cout << "Long Swap  ? " << ((ms_bSwapLong)?("TRUE"):("FALSE")) << "\n";
        cout << "Float Swap ? " << ((ms_bSwapFloat)?("TRUE"):("FALSE")) << "\n";
        cout << "Vax Float  ? " << ((ms_bVax2IEEE)?("TRUE"):("FALSE")) << "\n" << flush;
    }


    // If we are swapping longs,then we have other fields in the reflection list that 
    // need this treatment as well!

    if (ms_bSwapLong) {
        nNumRefs = nSwapLong(nNumRefs);
        nMaxString = nSwapLong(nMaxString);
        nExtraSize = nSwapLong(nExtraSize);
        nVersion = nSwapLong(nVersion);
        for (nx=0;nx<3;nx++) 
            a3nFields[nx] = nSwapLong(a3nFields[nx]);
    }

    // Check the integrity.

    if (nStat) {
        int nSizeFloat = (int) cSizeFloat;
        int nSizeInt = (int) cSizeInt;

        if ((nSizeFloat!=4) && (nSizeFloat!=8))
            nStat++;
        if (nSizeInt!=4)
            nStat++;
        if ((nMaxString<1) && (a3nFields[2]))
            nStat++;
        
        if (nSizeInt != sizeof(int)) {
            cout << "ERROR: Binary reflection list is incompatible with integers\n" << flush;
            nStat++;
        }
        if (nNumRefs < 0) {
            cout << "ERROR: Bad number of reflections in binary list\n" << flush;
            nStat++;
        }
        if (nExtraSize < 0) {
            cout << "ERROR: Bad 'extra size' field\n" << flush;
            nStat++;
        }
    }

    if( !nStat )
    {
        Cstring sTemp("");

        // Read any data that follows the header, yet preceeds the reflection data.
        if( nExtraSize > 0 && ppcExtra )
        {
            if( nBufferSizeIn - nCharsRead - 1 < nExtraSize ) 
                // "- 1", because the last character in the buffer is 0 and did not come from the file
            {
                int     nNewBufferSize = nCharsRead + nExtraSize + 1; // one extra character to have 0 at the end just in case

                // Re-allocate buffer
                delete pBufferIn;
                *ppcBuffer = new char[nNewBufferSize];
                memset(*ppcBuffer, 0, nNewBufferSize);
                
                // Re-read file
                FILE*   pFIn = fopen(sFilename.string(),"rb");
                if( NULL == pFIn )
                {
                    cout << "ERROR: could not open " << sFilename << '\n';
					cout << strerror(errno) << std::endl; //mrp
                    return 1;
                }
    
                if( 1 != fread(*ppcBuffer, nNewBufferSize - 1, 1, pFIn) )
                    nStat = 1;
                fclose(pFIn);
            }

            *ppcExtra = new char[nExtraSize];
            memset(*ppcExtra, 0, nExtraSize);
            memcpy(*ppcExtra, &((*ppcBuffer)[nCharsRead]), nExtraSize);
            nCharsRead += nExtraSize;
        }             
    }

    return nStat;
}


int tagBinary_Ref_Header::nWriteHeader(FILE* pFOut) {
    fprintf(pFOut,"\n\n%s",ms_pcId);
    fwrite((char*) this,sizeof(tagBinary_Ref_Header),1,pFOut);
    return 0;
}


int tagBinary_Ref_Header::nLoadHeader(Creflnlist& oReflnlist,const char* pcSelectIn) {
    int nx;
    int nLargest;
    int nRef;
    bool bWrite;
    char* pcComp;   // Pointer to the selectIn.

    nNumRefs= oReflnlist.nGetNumReflns();
    cSizeFloat = (char) sizeof(CReffloat);
    cSizeInt = (char) sizeof(int);

    // Search through all of the strings to find the longest one.  
    nLargest=0;
    for (nx=0;nx<oReflnlist.m_nCstringReflnFields;nx++) {
        for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) {
            nLargest = max((int) (oReflnlist[nRef].sGetField(nx).length()+1),nLargest);
        }
    }
    nMaxString = max(100,nLargest);



    a3nFields[0] = 0;
    a3nFields[1] = 0;
    a3nFields[2] = 0;

    pcComp = (char*) pcSelectIn;

    ms_sFields="";

    // Find the number and types of fields that we are writing.
    for (nx = 0; nx < 3; nx++) {
        bWrite = (pcComp==NULL);
        if (!bWrite) 
            bWrite = (*pcComp++ == 0);
        if (bWrite) {
            a3nFields[0]++;
            ms_sFields += oReflnlist.sGetFieldName(nx,eReflnField_int_type);
            ms_sFields += " ";
        }
    }

    for (nx = 3; nx < oReflnlist.m_nIntReflnFields; nx++) {
        bWrite = (pcComp==NULL);
        if (!bWrite) 
            bWrite = (*pcComp++ == 0);
        if (bWrite) {
            a3nFields[0]++;
            ms_sFields += oReflnlist.sGetFieldName(nx,eReflnField_int_type);
            ms_sFields += " ";
        }
    }
    for (nx = 0; nx < oReflnlist.m_nFloatReflnFields; nx++) {
        bWrite = (pcComp==NULL);
        if (!bWrite) 
            bWrite = (*pcComp++ == 0);
        if (bWrite) {
            a3nFields[1]++;
            ms_sFields += oReflnlist.sGetFieldName(nx,eReflnField_float_type);
            ms_sFields += " ";
        }
    }
    for (nx = 0; nx < oReflnlist.m_nCstringReflnFields; nx++) {
        bWrite = (pcComp==NULL);
        if (!bWrite) 
            bWrite = (*pcComp++ == 0);
        if (bWrite) {
            a3nFields[2]++;
            ms_sFields += oReflnlist.sGetFieldName(nx,eReflnField_Cstring_type);
            ms_sFields += " ";
        }
    }

    for (nx=0;nx<16;nx++)
    {
        _unionbrace._structbrace.a16cSampleFloat[nx]=(char) 0;
        _unionbrace._structbrace.a16cSampleInt[nx] = (char) 0;
    }
    
    *((CReffloat*) (&_unionbrace._structbrace.a16cSampleFloat[0])) = ms_fCodeFloat;
    *((int*) (&_unionbrace._structbrace.a16cSampleInt[0])) = ms_nCodeInt;

    nVersion = 0;

   return 0;
}

int tagBinary_Ref_Header::nWrite(FILE* pFOut,Crefln& oRefln,const char* pcSelectIn) {
    
    static char* pcBuf = NULL;
    static int  nBufSize = 0;
    char* pcPoint;
    int nx,ny;
    int nTotalStringBytesWritten = 0;
    int nSizeRecord;
    
    Creflnlist& oReflnlist = *(oRefln.m_poReflnlist);
    
    nSizeRecord = a3nFields[0]*((int) cSizeInt) + a3nFields[1]*((int) cSizeFloat) + a3nFields[2]*nMaxString;
    
    if ((pcBuf==NULL) || (nBufSize<nSizeRecord)) {
        if (pcBuf)
            delete[] ((int*) pcBuf);
        pcBuf = (char*) new int[nSizeRecord/sizeof(int) + 1];
        nBufSize=nSizeRecord;
    }
    
    pcPoint = pcBuf;
    // Write out all fields
    
    if (pcSelectIn == NULL) {
        for (nx = 0; nx < 3; nx++) {                    
            *((int*) pcPoint) = oRefln.nGetField(nx);
            pcPoint = (char*) (((int*) pcPoint)+1);
            
        }
        for (nx = 3; nx < oReflnlist.m_nIntReflnFields; nx++) {
            *((int*) pcPoint) = oRefln.nGetField(nx);
            pcPoint = (char*) (((int*) pcPoint)+1);
            
        }
        for (nx = 0; nx < oReflnlist.m_nFloatReflnFields; nx++) {
            *((CReffloat*) pcPoint) = oRefln.fGetField(nx);
            pcPoint = (char*) (((CReffloat*) pcPoint)+1);
            
        }
        for (nx = 0 , nTotalStringBytesWritten=0; nx < oReflnlist.m_nCstringReflnFields; nx++) {
            for (ny = 0; ny< oRefln.sGetField(nx).length(); ny++) {
                *pcPoint = oRefln.sGetField(nx).string()[ny];
                pcPoint++;
                nTotalStringBytesWritten++;
            }
            *pcPoint = (char) 0;
            pcPoint++;
            nTotalStringBytesWritten++;
        }
    } else if (oRefln.m_nSelect == pcSelectIn[oReflnlist.m_nTotalFieldsPlus1-1])
    {      

        const char *pcIn = pcSelectIn;        
        for (nx = 0; nx < 3; nx++) {                    
            if (0==*pcIn++) {
                *((int*) pcPoint) = oRefln.nGetField(nx);
                pcPoint = (char*) (((int*) pcPoint)+1);
            }
        }
        for (nx = 3; nx < oReflnlist.m_nIntReflnFields; nx++) {           
            if (0==*pcIn++) {
                *((int*) pcPoint) = oRefln.nGetField(nx);
                pcPoint = (char*) (((int*) pcPoint)+1);
            }           
        }
        for (nx = 0; nx < oReflnlist.m_nFloatReflnFields; nx++) {
            if (0==*pcIn++) {
                *((CReffloat*) pcPoint) = oRefln.fGetField(nx);
                pcPoint = (char*) (((CReffloat*) pcPoint)+1);
            }
        }
        for (nx = 0, nTotalStringBytesWritten = 0; nx < oReflnlist.m_nCstringReflnFields; nx++) {
            if (0==*pcIn++) {
                for (ny = 0; ny< oRefln.sGetField(nx).length(); ny++) {
                    *pcPoint = oRefln.sGetField(nx).string()[ny];
                    pcPoint++;
                    nTotalStringBytesWritten++;
                }
                *pcPoint = (char) 0;
                pcPoint++;
                nTotalStringBytesWritten++;
            }
        }
    }
    // Some computers require that the data be on a 4 byte boundry. 
    // We pad the strings so that this is so.
    nTotalStringBytesWritten = ((3*nTotalStringBytesWritten) % 4);
    for (nx = 0; nx < nTotalStringBytesWritten; nx ++) {
        *pcPoint = (char) 0;
        pcPoint++;
    }  
    
    fwrite(pcBuf,pcPoint-pcBuf,1,pFOut);
    return 0;
}


int tagBinary_Ref_Header::nRead(char* pcBuf,Crefln& oRefln) {
    int nx,ny,nz;
    int nStat = 0;

  static Cstring sTemp;
  char*  pcBufIn;

  pcBufIn = pcBuf;

  Creflnlist& oReflnlist = *(oRefln.m_poReflnlist);
  
  // Read all the integer fields
  for (nx = 0; (nx < oReflnlist.m_nIntReflnFields) && (0 == nStat); nx++) 
    {
      int nTemp;

      nTemp = *((int*) pcBuf);
      pcBuf = (char*) (((int*) pcBuf)+1);
      if (ms_bSwapLong) 
          nTemp = nSwapLong(nTemp);
      oRefln.vSetField(nx, nTemp);
    }

  // Read all the float fields

  for (nx = 0; (nx < oReflnlist.m_nFloatReflnFields) && (0 == nStat); nx++) {
      float fTemp;
      double fDoubleTemp;

      if (((int) cSizeFloat) == 8) {
          fDoubleTemp = *((double*) pcBuf);
          pcBuf = (char*) (((double*) pcBuf)+1);
          oRefln.vSetField(nx, (float) fDoubleTemp);
      } else {
          fTemp = *((float*) pcBuf);
          pcBuf = (char*) (((float*) pcBuf)+1);
          if (ms_bVax2IEEE) 
              fTemp = fVAXtoIEEE(&fTemp);
          if (ms_bSwapFloat)
              fTemp = fSwapFloat(fTemp);          
          oRefln.vSetField(nx, fTemp);
      }
  }

  // Read all the Cstring fields


  for (nx = 0, nz=0; (nx < oReflnlist.m_nCstringReflnFields) && (0 == nStat); nx++) {
      sTemp = pcBuf;
      ny = sTemp.length()+1;
      nz += ny;
      pcBuf += ny;
      oRefln.vSetField(nx,sTemp);
  }
  // skip any padding.
  pcBuf += ((3*nz) % 4);

  return (pcBuf - pcBufIn);
}


void Creflnlist::vWritePartial(const Cstring & sFilename)
{
  //Crefln *poRefln;
  //int i;
  //int nStat;
  int nFirstLocal;
  int nLastLocal;
  //bool bNewIndex = FALSE;

  //FILE* pFOut = fopen(sTransSymbol(sFilename).string(),"wb");
  nFirstLocal = 0;
  nLastLocal  = (int)((float)m_nNumReflns*m_fPercentRefToDisplay);
  int* iSortIndex = pnGetSortIndex();

  nWrite(sFilename, iSortIndex, NULL, NULL, nFirstLocal, nLastLocal);

  /*nStat = 0;
  if (pFOut)
    {
      nWriteHeader(pFOut, NULL,false);
      for (i = nFirstLocal; (i < nLastLocal) && (0 == nStat); i++)
        {
          int idx = i;
          if(iSortIndex) {
             idx = iSortIndex[i];
          }
          poRefln = m_ppoTheReflns[idx];
          poRefln->vWrite(pFOut, NULL);
          if (ferror(pFOut))
            {
              fclose(pFOut);
              pFOut = NULL;
              nStat = 1;
            }
        }
      if (pFOut)
        fclose(pFOut);
      m_bNeedsSaving = 0;
    }*/
} 
int Creflnlist::nUpdateCrystal(Cstring& sFileName, Ccrystal* poCrystal, Cstring* psCommentLog)
{
    FILE* pFIn;
    FILE* pFOut;
    int nx;
    char pcBuffer[200];

    vDeleteAll();
    if (nRead(sFileName,0)) {
        return 1;
    }
    if (poCrystal)
        nPutCrystal(*poCrystal);
    if (psCommentLog) 
        nLogComment(*psCommentLog,FALSE);

    pFIn = fopen(sFileName.string(),"rb");
    if (!pFIn) 
        return 1;

    if (m_eWriteBinary == eReflnoutput_text) {
        for (nx=0;nx<m_nBytesOrLinesReadInHeader;nx++) {
            if ((!fgets(pcBuffer,200-1,pFIn)) || (pcBuffer[strlen(pcBuffer)-1]!='\n')) {
                fclose(pFIn);
                return 1;
            }
        }
    } else {        
        for (nx=0;nx<m_nBytesOrLinesReadInHeader;nx++) 
            fread(&pcBuffer[0],1,1,pFIn);
            
    }
    
    // Find how much data there remains to be read.
    int nPosInit;
    int nPosEnd;
    int nSize;
    int nMemorySize;
    char* pcBuf;
    
    
    nPosInit = ftell(pFIn);
    fseek(pFIn,0,SEEK_END);
    nPosEnd = ftell(pFIn);
    fseek(pFIn,nPosInit,SEEK_SET);
    nMemorySize = nSize = nPosEnd - nPosInit;
    // read in the data as a block.       
    pcBuf = new char[nMemorySize]; 
    fread(pcBuf,nSize,1,pFIn);
    fclose(pFIn);
    
    // Open the old file name for writing.
    pFOut = fopen(sFileName.string(),"wb");
    if (!pFOut) {
        delete[] pcBuf;
        return 1;
    }
       
    nWrite(sFileName, NULL, NULL, pFOut);
    fwrite(pcBuf,nSize,1,pFOut);
    fclose(pFOut);    
    delete[] pcBuf;
    
    return 0;
}


int
Creflnlist::nConvertFromAnom(void)
{
  // Take a list with anomalous I+ and I- fields and convert it those fields
  // to hkl, Int, sigInt
  // -h-k-l, Int, sigInt

  // Problem, what about centric reflections!

  int     i, j;             // Loop counter
  Crefln  oRefln(this);     // A reflection object for filling and inserting

  // Look for necessary fields

  if (   (0 > m_nFI_fIntensityPlus)
      || (0 > m_nFI_fIntensityMinus)
      || (0 > m_nFI_fSigmaIPlus)
      || (0 > m_nFI_fSigmaIMinus) )
    {
      cout << "ERROR in Creflnlist::nConvertFromAnom(), required fields missing.\n" << endl;
      return (-1);
    }

  Crefln **ppoNewReflns;
  int    nNumAlloc;
  nNumAlloc    = 2 * nGetNumReflns();
  ppoNewReflns = new Crefln* [nNumAlloc];
  for (i = 0; i < nNumAlloc; i++)
    *(ppoNewReflns+i) = NULL;

  j = 0;
  Crefln *poRefln;
  int nH, nK, nL;
  for (i = 0; i < nGetNumReflns(); i++)
    {
      poRefln = *(m_ppoTheReflns+i);
      nH = poRefln->nGetH();
      nK = poRefln->nGetK();
      nL = poRefln->nGetL();
      ppoNewReflns[j] = new Crefln(this, nH, nK, nL);
      if (   (-1.0 != poRefln->fGetField(m_nFI_fIntensityPlus))
          && (-1.0 != poRefln->fGetField(m_nFI_fSigmaIPlus)) )
        {
          ppoNewReflns[j] = new Crefln(this, nH, nK, nL);
          ppoNewReflns[j]->vSetIntensity(poRefln->fGetField(m_nFI_fIntensityPlus));
          ppoNewReflns[j]->vSetSigmaI(poRefln->fGetField(m_nFI_fSigmaIPlus));

          j++;
        }
      if (   (-1.0 != poRefln->fGetField(m_nFI_fIntensityMinus))
          && (-1.0 != poRefln->fGetField(m_nFI_fSigmaIMinus)) )
        {
          ppoNewReflns[j] = new Crefln(this, -nH, -nK, -nL);
          ppoNewReflns[j]->vSetIntensity(poRefln->fGetField(m_nFI_fIntensityMinus));
          ppoNewReflns[j]->vSetSigmaI(poRefln->fGetField(m_nFI_fSigmaIMinus));

          j++;
        }
    }

  for (i = 0; i < m_nNumReflns; i++)
    {  // Delete each reflection
      if (NULL != m_ppoTheReflns[i])
        {
          delete m_ppoTheReflns[i];
          m_ppoTheReflns[i] = NULL;
        }
    }
  if (0 < m_nNumAlloc)
    {
      delete [] m_ppoTheReflns;
      m_ppoTheReflns = NULL;
    }
  m_ppoTheReflns = ppoNewReflns;
  m_nNumReflns   = j;
  m_nNumAlloc    = nNumAlloc;
  vUpdateFieldIndices();

/*
NOT YET IMPLEMENTED
  // Deselect the 4 float fields for possible no output later

  for (j = 0; j < m_nFloatReflnFields; j++)
    {
      if (      (ms_sfIntensityPlus  == *(m_psFloatNames+j))
                || (ms_fSigmaIPlus      == *(m_psFloatNames+j))
                || (ms_sfIntensityMinus == *(m_psFloatNames+j))
                || (ms_sfSigmaIMinus    == *(m_psFloatNames+j))      )
        {
          // set
        }
    }
*/
  return (0);
}
int Creflnlist::nList(Cimage_header& oHeader)
{
    double f0;
    int nComponents;
    int nTwinLaw;
    int nComponentsFromCrystal;
    int nStat;
    Ccrystal* m_poCrystal;

    printf("\n");
    printf("Reflection listing\n");
    printf("---------------------------------------------------------------------------\n");
    printf(" Component Num Indexed   <MM-RMS>  Spgp/Latt         a         b         c \n");
    printf("  Rotation <Intensity>  <DEG-RMS>     Volume     alpha      beta     gamma \n");
    printf("    Vector     <I/sig>  <HKL-RMS>  Mosaicity      rot1      rot2      rot3 \n");
    printf("---------------------------------------------------------------------------\n");

    if (nCalcRecipLatticePoints(oHeader,f0 = 10.0))
        return 1;
    nComponents = 0;
    nStat = 0;
    do {
        m_poCrystal = new Ccrystal(oHeader,nComponents+1); 
        if ((!m_poCrystal->bIsAvailable()) && (nComponents==0))
            break;
        else if (!m_poCrystal->bIsAvailable())
            break;
        else {
            for (nTwinLaw = 0; nTwinLaw <= m_poCrystal->m_nTwinLaws; nTwinLaw++) {
                double afValues[g_nResidInfoRMSEntries];
                double fRot,fRotSig;
                double a3fRotVector[3];
                int    a3nRot[3];
                double a3x3fOrient[3][3];
                double a3x3fTMat[3][3];
                double a3x3fTOrient[3][3];
                double a3fOrient[3];
                char pcTwinLawBuf[50];


                Ccrystal::nGetSetResidValues(oHeader,nComponents + 1,nTwinLaw,&afValues[0],0);
                m_poCrystal->nCalcOrientMatrix(1.0);
                m_poCrystal->nGetTwinLawRotAndVector(nTwinLaw,3,&fRot,&fRotSig,&a3fRotVector[3],&a3nRot[0]);
                m_poCrystal->nCalcRotMatrix();
                m_poCrystal->vGetRotMatrix(&a3x3fOrient[0][0]);
                m_poCrystal->nGetTwinLaw(nTwinLaw,&a3x3fTMat[0][0],true);
                vMulMat3DMat3D(a3x3fTMat,a3x3fOrient,a3x3fTOrient);
                vDeconvMat3D3XYZ(a3x3fTOrient,&a3fOrient[0],&a3fOrient[1],&a3fOrient[2]);


                if ((nTwinLaw) || (nComponents))
                    printf("\n");

                printf(" %2d law %2d %11d %10.2lf    %4d %c%c %9.4lf %9.4lf %9.4lf\n",
                    nComponents + 1,nTwinLaw + 1,(int) afValues[g_nResidInfoNumIndexed],afValues[g_nResidInfoRMSMm],
                    m_poCrystal->m_poSpacegroup->nGet(),m_poCrystal->m_poSpacegroup->cGetClass(),g_cpSpaceGroupNames[m_poCrystal->m_poSpacegroup->nGet()-1][0],
                    m_poCrystal->fGetCell(0),m_poCrystal->fGetCell(1),m_poCrystal->fGetCell(2));
                if (nTwinLaw == 0)
                    pcTwinLawBuf[0] = 0;
                else
                    sprintf(pcTwinLawBuf,"%9.2lf",fRot);
                printf(" %9s %11.0lf %10.2lf  %9.0lf %9.4lf %9.4lf %9.4lf\n",
                    pcTwinLawBuf,afValues[g_nResidInfoIntensity],afValues[g_nResidInfoRMSRot],m_poCrystal->fCalcVolume(),
                    m_poCrystal->fGetCell(3),m_poCrystal->fGetCell(4),m_poCrystal->fGetCell(5));
                if (nTwinLaw == 0)
                    pcTwinLawBuf[0] = 0;
                else
                    sprintf(pcTwinLawBuf,"[%2d%2d%2d]",a3nRot[0],a3nRot[1],a3nRot[2]);
                printf(" %9s %11.2lf %10.2lf  %9.1lf %9.1lf %9.1lf %9.1lf\n",
                    pcTwinLawBuf,afValues[g_nResidInfoIoverSig],afValues[g_nResidInfoRMSHKL],
                    m_poCrystal->fGetMosaicity(),a3fOrient[0],a3fOrient[1],a3fOrient[2]);
            }
            nComponentsFromCrystal = m_poCrystal->m_nTwins;
        }
        delete m_poCrystal;
        m_poCrystal = NULL;
        nComponents++;
    } while ((!nStat) && (nComponents<nComponentsFromCrystal));

    printf("---------------------------------------------------------------------------\n");
    
    if (m_poCrystal)
        delete m_poCrystal;
    m_poCrystal = NULL;


    return 0;
}


int Creflnlist::nCalcRecipLatticePoints(Cimage_header& oHeader,double& fUserInput)
{
  int   i, nStat, nDetNum, nDetNumLast, nRef;
  float  a3fFloatBuf[3];
  float  f1MM, f2MM, f3MM;
  float  fCalcPx0, fCalcPx1;
  double fObsPx0,fObsPx1;
  double fDD[3][3],fDDInv[3][3],fDN[3];         // Detector ....
  double fVDCO[3];                  // Virtual detector observed coordinates
  double fTS[3];
  double fLorentz;                  // Mosaicity variable.
  double fDStarSq;                  // Mosaicity variable.
  double fDSCOTH;                   // Mosaicity variable.
  double a3fTempVec[3];

  double fRecipShiftRot[3];         // Recip. shift vector rotated
  double fXR[3];                    // Recip. lattice coords in diffracting condition
  double fX[3];                     // Unrotated Recip latt coords
  double fXRD[3];                   // Recip. lattice derivative coords in diffracting condition.
  double fXD[3];                    // Unrotated Recip. lattice derivative coords.
  double fObsRot;                   // Obs Rotation angle
  double fCalcRot;                  // Calc Rotation angle
  double fRotVec[3];
  double fInvRot[3][3];
  double fGonMatrix[3][3];
  double fRotGonMatrix[3][3];
  double fInvGonMatrix[3][3];
  double fCrysOrient[20][g_nMaxTwinLaws + 1][3][3];         // Orient matrix with twin law applied.
  double fCrysOrientInv[20][g_nMaxTwinLaws + 1][3][3];      // Inverse of above.
  double fRecipShift[20][g_nMaxTwinLaws + 1][3];            // Reciprocal shift for each twin law
  double fMosaicity[20][3];                                 // Mosaicity for each component.
  int    anNumTwinLaws[20];                                 // # of twin laws for each component.
  double afRMSMM[20][g_nMaxTwinLaws + 1];                   // These are used for storing information on indexing.
  double afRMSDEG[20][g_nMaxTwinLaws + 1];
  double afRMSHKL[20][g_nMaxTwinLaws + 1];
  double afIntensity[20][g_nMaxTwinLaws + 1];
  double afIoverSig[20][g_nMaxTwinLaws + 1];
  int    anIntensityContrib[20][g_nMaxTwinLaws + 1];
  int    anContrib[20][g_nMaxTwinLaws + 1];
  int    nTotalContrib;



  double fDet;
  double a3fS0[3];
  double fWavelength;
  double fMinTwinIDTolerance;       // Nominally, the user input.
  double fRejectFraction = 0.0;     // The user input if it is <0.
  int   nFI0, nFI1, nFI2,nFI3;
  int   nComp,nTwinLaw;
  int   nComponents;
  int   nComponentsFromCrystal;
  bool  bComputeTwinID;
  bool  bEllipsoidIntersects;
  bool  bEllipsoidInfoAvailable;

  Crefln *poRefln;                 // Convenience pointer

  if (fUserInput<=-1.0) {
      bComputeTwinID = false;
  } else if (fUserInput>=0.0) {
      bComputeTwinID = true;
      fRejectFraction = 0.0;
      fMinTwinIDTolerance = fUserInput;
  } else {
      fRejectFraction = -fUserInput;
      fMinTwinIDTolerance = 100.0;
      bComputeTwinID = true;
  }

  if (0 > m_nFI_fObsPx0)
    {
      nFI0 = m_nFI_fCalcPx0;
    }
  else
    {
      nFI0 = m_nFI_fObsPx0;
    }
  if (0 > m_nFI_fObsPx1)
    {
      nFI1 = m_nFI_fCalcPx1;
    }
  else
    {
      nFI1 = m_nFI_fObsPx1;
    }
  if (0 > m_nFI_fObsRotMid)
    {
      nFI2 = m_nFI_fCalcRotMid;
    }
  else
    {
      nFI2 = m_nFI_fObsRotMid;
    }
  if (0 > m_nFI_fObsRotWidth)
  {
      nFI3 = m_nFI_fCalcRotWidth;
  }
  else
  {
      nFI3 = m_nFI_fObsRotWidth;      
  }

  if ( (0 > nFI0) || (0 > nFI1) || (0 > nFI2) )
    {
      printf("ERROR, missing spot centroids from reflection list!\n");
      return (1);
    }

  // Expand to get fields that we need.
  nExpandGetField(ms_sfRecipCoord0);
  nExpandGetField(ms_sfRecipCoord1);
  nExpandGetField(ms_sfRecipCoord2);

  if (bComputeTwinID) 
      nExpandGetField(ms_snTwinID);

  // Load the various objects from the header.
  Ccrystal* m_poCrystal             = NULL;
  Cgoniometer* m_poCrysGonio        = NULL;
  Crotation* m_poRotation           = NULL;
  Csource* m_poSource               = NULL;
  Cdetector** m_ppoDetector         = NULL;
  Cstring* m_psDetectorNames        = NULL;
  int      m_nNumDetectors          = 0;
  float*   pfResidual               = NULL;

  m_poCrysGonio = new Cgoniometer(oHeader, Ccrystal::ms_sCrystalPrefix);
  m_poRotation    = new Crotation(oHeader);
  m_poSource = new Csource(oHeader);

  if (fRejectFraction>0.0)
      pfResidual = new float[nGetNumReflns()];
  
  nStat = oHeader.nGetValue(Cdetector::ms_sDetectorNumber, &m_nNumDetectors);  
  if (0 == nStat) {
    m_psDetectorNames = new Cstring [m_nNumDetectors];
    nStat = oHeader.nGetValue(Cdetector::ms_sDetectorNames,
      m_nNumDetectors, m_psDetectorNames);
  }
  if (0 == nStat)
    {
      m_ppoDetector = new Cdetector* [m_nNumDetectors];
      for (i = 0; i < m_nNumDetectors; i++) 
        {
          m_ppoDetector[i] = new Cdetector (oHeader,
                                            m_psDetectorNames[i],TRUE, FALSE);
          if (!m_ppoDetector[i]->bIsAvailable())
            {
              printf("Could not load detector.\n");
              nStat = 1;
              break;
            }
        }
    }  
  
  // Load the crystal orientation matrices.
  if (bComputeTwinID) {
      nComponents = 0;
      nStat = 0;
      do {
          m_poCrystal = new Ccrystal(oHeader,nComponents+1); 
          if ((!m_poCrystal->bIsAvailable()) && (nComponents==0))
              break;
          else if (!m_poCrystal->bIsAvailable())
              break;
          else {
              for (nTwinLaw = 0; nTwinLaw <= m_poCrystal->m_nTwinLaws; nTwinLaw++) {
                  double a3x3fTemp[3][3];
                  double a3x3fTemp2[3][3];
                  double a3x3fTwinLaw[3][3];

                  m_poCrystal->nGetTwinLaw(nTwinLaw,&a3x3fTwinLaw[0][0],true);
                  m_poCrystal->nCalcOrientMatrix(m_poSource->fGetWavelength()); 
                  m_poCrystal->vGetOrientMatrix(& a3x3fTemp[0][0]);
                  vMulMat3DMat3D(a3x3fTwinLaw,a3x3fTemp,a3x3fTemp2);
                  vCopyMat3D(&a3x3fTemp2[0][0],&fCrysOrient[nComponents][nTwinLaw][0][0]);
                  fInvMat3D(&fCrysOrient[nComponents][nTwinLaw][0][0],&fCrysOrientInv[nComponents][nTwinLaw][0][0]);
                  m_poCrystal->nGetRecipShift(nTwinLaw,&a3x3fTemp[0][0]);
                  vCopyVec3D(&a3x3fTemp[0][0],&fRecipShift[nComponents][nTwinLaw][0]);
                  m_poCrystal->vGetMosaicity(&fMosaicity[nComponents][0]);
                  
                  afRMSMM[nComponents][nTwinLaw] = 0.0;
                  afRMSDEG[nComponents][nTwinLaw] = 0.0;
                  afRMSHKL[nComponents][nTwinLaw] = 0.0;
                  afIntensity[nComponents][nTwinLaw] = 0.0;
                  afIoverSig[nComponents][nTwinLaw] = 0.0;
                  anIntensityContrib[nComponents][nTwinLaw] = 0;
                  anContrib[nComponents][nTwinLaw] = 0;
              }
              anNumTwinLaws[nComponents] = m_poCrystal->m_nTwinLaws;
              nComponentsFromCrystal = m_poCrystal->m_nTwins;
          }
          delete m_poCrystal;
          m_poCrystal = NULL;
          nComponents++;
      } while ((!nStat) && (nComponents<nComponentsFromCrystal));
      if (m_poCrystal)
          delete m_poCrystal;
      m_poCrystal = NULL;
      nTotalContrib = 0;
  }
  
  do {
      if ((!nStat) && (m_poCrysGonio->bIsAvailable()) && (m_poRotation->bIsAvailable()) &&
          (m_poSource->bIsAvailable())) {
          
          // Now calculate the reciprocal lattice coordinates of the reflections
          // from the detector coordinates!
          
          // Note:  In d*TREK coordinate system:
          //    S0 points from origin of reciprocal lattice towards origin of crystal
          //
          // xr = R G U B h
          //         -1   -1
          // xu =   G  * R  * (S + S0)  , unrotated reciprocal lattice coordinate
          // --     =    =     -   --
          //
          //  tS = DD * VDCO
          //   -   ==   ----
          //
          //  S = tS / |tS|,  assuming the Ewald sphere has been normalized to a radius
          //  -    -     -    of 1.
          //
          //  XR = S + S0
          //  --   -   --
          //
          //   X = InvG * InvR * XR     : rlp coords referenced to goniometer at 0,0,0
          //   -   ====   ====   --
          
          
          // Get crystal goniometer rotation matrix at datum
          
          m_poCrysGonio->vCalcGetRotMatrix(&fGonMatrix[0][0],m_poCrysGonio->nGetNum(m_poRotation->sGetName()));
          fDet = fInvMat3D(&fGonMatrix[0][0],&fInvGonMatrix[0][0]);
          if (m_poCrysGonio->nGetRotVector(m_poRotation->sGetName(),&a3fFloatBuf[0])) {
              cout << "ERROR, check that ROTATION_AXIS_NAME value matches a CRYSTAL_GONIO_NAMES value!\n";
              nStat =1;
              continue;
          }
          vCopyVec3D(&a3fFloatBuf[0],&fRotVec[0]);
          m_poRotation->vSetVector(&a3fFloatBuf[0]);
          
          // Get source vector, normalized
          
          m_poSource->vCalcGetS0(a3fS0);
          fWavelength = m_poSource->fGetWavelength();
          
          // Get Detector matrix and vector
          
          nDetNumLast = 0;
          m_ppoDetector[0]->nCalcGetDDDN(&fDD[0][0], &fDN[0]);
          m_ppoDetector[0]->vUpdateDDDN();
          fInvMat3D(&fDD[0][0], &fDDInv[0][0]);
          
          
          nRef  = nGetNumReflns();
          for (i = 0; i < nRef; i++) {

              // In this loop over reflections, we should probably reject
              //  reflections outside resolution and I/sigma cutoffs, too.
              
              // Get detector number for this reflection
              
              poRefln = poGetRefln(i);

              if (m_nFI_nDetNum>=0)
                  nDetNum = poRefln->nGetField(m_nFI_nDetNum);
              else
                  nDetNum = 0;
              if (nDetNum != nDetNumLast)
              {
                  m_ppoDetector[nDetNum]->nCalcGetDDDN(&fDD[0][0], &fDN[0]);
                  m_ppoDetector[nDetNum]->vUpdateDDDN();
                  nDetNumLast = nDetNum;
              }

              // See if the goniometer matrix needs to be recalculated.
              // This will only be the case if we have multiple scans with different
              // goniometer axes incoded by the m_nFI_fGonio?  fields of the reflection list.
              if (m_nFI_fGonio1>=0) {
                  int nGonioAxis;
                  bool bChanged;

                  if (m_nFI_nGonioRotAxis<0) {
                      printf("ERROR:  Goniometer rotation axis field is missing.\n");
                      nStat = 1;
                      continue;
                  }
                  
                  bChanged = false;
                  nGonioAxis =0;
                  if ((m_nFI_fGonio1>=0) && (ABS(poRefln->fGetField(m_nFI_fGonio1) - m_poCrysGonio->m_pfDatumValue[nGonioAxis])>0.01)) {
                      m_poCrysGonio->m_pfDatumValue[nGonioAxis] = poRefln->fGetField(m_nFI_fGonio1);
                      bChanged = true;
                  }
                  nGonioAxis++;
                  if ((m_nFI_fGonio2>=0) && (ABS(poRefln->fGetField(m_nFI_fGonio2) - m_poCrysGonio->m_pfDatumValue[nGonioAxis])>0.01)) {
                      m_poCrysGonio->m_pfDatumValue[nGonioAxis] = poRefln->fGetField(m_nFI_fGonio2);
                      bChanged = true;
                  }
                  nGonioAxis++;
                  if ((m_nFI_fGonio3>=0) && (ABS(poRefln->fGetField(m_nFI_fGonio3) - m_poCrysGonio->m_pfDatumValue[nGonioAxis])>0.01)) {
                      m_poCrysGonio->m_pfDatumValue[nGonioAxis] = poRefln->fGetField(m_nFI_fGonio3);
                      bChanged = true;
                  }
                  nGonioAxis++;
                  if ((m_nFI_fGonio4>=0) && (ABS(poRefln->fGetField(m_nFI_fGonio4) - m_poCrysGonio->m_pfDatumValue[nGonioAxis])>0.01)) {
                      m_poCrysGonio->m_pfDatumValue[nGonioAxis] = poRefln->fGetField(m_nFI_fGonio4);
                      bChanged = true;
                  }
                  nGonioAxis++;
                  if ((m_nFI_fGonio5>=0) && (ABS(poRefln->fGetField(m_nFI_fGonio5) - m_poCrysGonio->m_pfDatumValue[nGonioAxis])>0.01)) {
                      m_poCrysGonio->m_pfDatumValue[nGonioAxis] = poRefln->fGetField(m_nFI_fGonio5);
                      bChanged = true;
                  }
                  nGonioAxis++;
                  
                  if (bChanged)
                    {
                      m_poCrysGonio->vCalcGetRotMatrix(&fGonMatrix[0][0],poRefln->nGetField(m_nFI_nGonioRotAxis)-1);
                      fDet = fInvMat3D(&fGonMatrix[0][0],&fInvGonMatrix[0][0]);
                      m_poCrysGonio->nGetRotVector(poRefln->nGetField(m_nFI_nGonioRotAxis)-1,&a3fFloatBuf[0]);
                      vCopyVec3D(&a3fFloatBuf[0],&fRotVec[0]);
                      m_poRotation->vSetVector(&a3fFloatBuf[0]);
                  }
              }
                    
              // Get spot centroid in pixel coordinates, degrees
              // Use calculated position if observed position not available
              
              fObsPx0 = poRefln->fGetField(nFI0);
              fObsPx1 = poRefln->fGetField(nFI1);
              fObsRot = poRefln->fGetField(nFI2);
              
              // Calculate observed millimeter coordinates
              
              nStat = m_ppoDetector[nDetNum]->m_poSpatial->nPixeltoMM(fObsPx0, fObsPx1,
                  &f1MM, &f2MM, &f3MM);
              
              if (0 != nStat) {                 
                  fX[0] = 0.0;
                  fX[1] = 0.0;
                  fX[2] = 0.0;
                  poRefln->vSetField(m_nFI_fRecipCoord0, (float) fX[0]);
                  poRefln->vSetField(m_nFI_fRecipCoord1, (float) fX[1]);
                  poRefln->vSetField(m_nFI_fRecipCoord2, (float) fX[2]);
                  if (bComputeTwinID) {                      
                      poRefln->vSetField(m_nFI_fHKLResid,(float) 100.0);
                      poRefln->vSetField(m_nFI_nTwinID,(int) 0);
                  }
                  nStat = 0;
              } else {
                  
                  // A single goniometer datum and rot. axis is used for all reflections
                  // fGonMatrix is calculated above outside of reflection loop
                  // Notice how by using a negative rotation that we get the inverse!
                  
                  vConvRotVec3DMat3D(-fObsRot, fRotVec, &fInvRot[0][0]);
                  
                  // Calculate reciprocal lattice coordinates at diffracting position
                  
                  fVDCO[0] = f1MM;
                  fVDCO[1] = f2MM;
                  fVDCO[2] = f3MM;
                  vMulMat3DVec3D(fDD, fVDCO, fTS);
                  (void) fNormVec3D(fTS);
                  vAddVec3DVec3D(a3fS0, fTS, fXR);


                  // Compute a "Lorentz factor" and reject this reflection if a problem                  
                  vCross3D(fRotVec, a3fS0, a3fTempVec);               // fTS is a temp var here
                  fLorentz = 1.0/max(ABS(fDot3D(a3fTempVec, fXR)),1e-10);                
                  fDStarSq = fDot3D(fXR,fXR);


                  // Obtain the derivative vector.
                  vMulVec3DScalar(fRotVec,fDot3D(fXR,fRotVec),fTS);
                  vSubVec3DVec3D(fXR,fTS,fTS);
                  vCross3D(fRotVec,fTS,fXRD);
                  
                  vMulMat3DVec3D(fInvRot, fXR, fTS);      // fTS is a temp var here
                  vMulMat3DVec3D(fInvGonMatrix, fTS, fX); // fTS is a temp var here

                  vMulMat3DVec3D(fInvRot, fXRD, fTS);      // fTS is a temp var here
                  vMulMat3DVec3D(fInvGonMatrix, fTS, fXD); // fTS is a temp var here                  

                  // Resolution is wavelength / |XR|
                  
                  fDet = fLenVec3D(fXR);
                                                                             
                  // Put fX in refln fields overwriting fXIn if it is there

                  poRefln->vSetField(m_nFI_fRecipCoord0, (float) fX[0]);
                  poRefln->vSetField(m_nFI_fRecipCoord1, (float) fX[1]);
                  poRefln->vSetField(m_nFI_fRecipCoord2, (float) fX[2]);

                  // Compute the Derivative of fX w.r.t. rotation.  This is used in indexing.
                  if ((m_nFI_fRecipCoordD0>=0) && (m_nFI_fRecipCoordD1>=0) && (m_nFI_fRecipCoordD2>=0)) {
                      // Store XD derivative vector.
                      poRefln->vSetField(m_nFI_fRecipCoordD0, (float) fXD[0]);
                      poRefln->vSetField(m_nFI_fRecipCoordD1, (float) fXD[1]);
                      poRefln->vSetField(m_nFI_fRecipCoordD2, (float) fXD[2]);
                  }
                  
                  if (bComputeTwinID) {
                      double a3x3fTempMat1[3][3];
                      double a3fTempVec1[3];
                      double a3x3fRotGonMatrix[3][3];
                      double a3x3fInvRotGonMatrix[3][3];
                      double a3fHKL[3];
                      int   a3nHKL[3];
                      int nBestTolerance;
                      int nBestTwinLaw;
                      double fBestHKLErr;
                      double fBestMMErr;
                      double fBestRotErr;
                      double fHKLErr,fMMErr,fRotErr;
                      double f0;
                      int nx;
                      vConvRotVec3DMat3D(fObsRot, fRotVec, & a3x3fTempMat1[0][0]);
                      vMulMat3DMat3D(a3x3fTempMat1, fGonMatrix, a3x3fRotGonMatrix);
                      vCopyMat3D(&a3x3fRotGonMatrix[0][0],&a3x3fInvRotGonMatrix[0][0]);
                      vTranMat3D(a3x3fInvRotGonMatrix);
                      
                      // Calculate the fractional hkl entries.
                      //              -1        -1
                      //       R = (TUB)  * ( PHR  * (VNTS + S0))
                      
                      fBestHKLErr = 10.0;
                      nBestTolerance = -1;
                      nBestTwinLaw = 0;
                      for (nComp=0;nComp<nComponents;nComp++) {
                          for (nTwinLaw = 0; nTwinLaw <= anNumTwinLaws[nComp]; nTwinLaw++) {
                              vMulMat3DVec3D(a3x3fInvRotGonMatrix, fXR, a3fTempVec1);
                              vMulMat3DVec3D(fCrysOrientInv[nComp][nTwinLaw], a3fTempVec1,a3fHKL);
                              fHKLErr = 0.0;
                              for (nx=0;nx<3;nx++) {
                                  if (a3fHKL[nx] - floor(a3fHKL[nx])<ceil(a3fHKL[nx])-a3fHKL[nx]) 
                                      a3nHKL[nx] = (int) floor(a3fHKL[nx]);
                                  else
                                      a3nHKL[nx] = (int) ceil(a3fHKL[nx]);
                                  
                                  f0 = a3fHKL[nx] - a3nHKL[nx];                              
                                  f0 *= f0;
                                  fHKLErr += f0;
                              }

                              bEllipsoidIntersects = false;
                              bEllipsoidInfoAvailable = false;
                              // If spot centroid information is available...

                              fRotErr = 0.0;
                              fMMErr = 0.0;

                              if ((m_nFI_fObsRotMid>=0) && 
                                  (m_nFI_fObsRotWidth>=0) && 
                                  (nFI3>=0)) {
                                  
                                  double a3fVTS[3];
                                  double a3fVDCC[3];
                                  double a3fXRC[3];                                  
                                  double a3fXRCR[3];
                                  
                                  // Then we project the spot back onto the plate...
                                  for (nx = 0; nx < 3; nx++) {
                                      a3fHKL[nx] = a3nHKL[nx];
                                  }    
                                  vMulMat3DVec3D(fCrysOrient[nComp][nTwinLaw],a3fHKL,a3fTempVec1);            
                                  vMulMat3DVec3D(a3x3fRotGonMatrix, a3fTempVec1, a3fXRC);
                                  
                                  fRotErr = fCalcRotationOffset(fRotVec,a3fS0,a3fXRC,a3fXRCR);
                                  if (fRotErr >= -999.0) {
                                      
                                      
                                      fCalcRot = fObsRot + fRotErr;
                                     
                                      // Calculate the rotated X-calc
                                      vConvRotVec3DMat3D(fRotErr,&fRotVec[0],&a3x3fTempMat1[0][0]);
                                      vMulMat3DVec3D(a3x3fTempMat1,a3fXRC,a3fXRCR);
                                      
                                      // This code transforms X-calc (rotated into difracting condition) into pixel coordinates
                                      // by first: adding in the S0 vector
                                      vSubVec3DVec3D(a3fXRCR, a3fS0, a3fTempVec1);
                                      // next:  accounting for the precessing reciprocal shift. (we recalculate a3fRecipVec)
                                      vMulMat3DVec3D(fRotGonMatrix,&fRecipShift[nComp][nTwinLaw][0],&fRecipShiftRot[0]);
                                      // lastly:  projecting onto the detector plate (taking geometry into account)
                                      if (!m_ppoDetector[0]->nScaleSOnPlateItr(a3fTempVec1,a3fVTS,fRecipShiftRot))
                                      {
                                          // Calculate the virtual detector coordinates by using the inverse DD matrix.
                                          
                                          vMulMat3DVec3D(fDDInv, a3fVTS, a3fVDCC);
                                          
                                          if (!m_ppoDetector[0]->m_poSpatial->nMMtoPixel(a3fVDCC[0],
                                              a3fVDCC[1],a3fVDCC[2],&fCalcPx0, &fCalcPx1)) {

                                              fMMErr = (f1MM - a3fVDCC[0])*(f1MM - a3fVDCC[0]) + (f2MM - a3fVDCC[1])*(f2MM - a3fVDCC[1]);
                                              
                                              double a2x2fFindEllipseA[2][2];
                                              double a2fFindEllipseb[2];
                                              double fFindEllipsec;
                                              double fFindRotMid;
                                              double fFindRotWidth;
                                              double fFindRotStart;
                                              double fFindRotEnd;
                                              double fCalcRotStart;
                                              double fCalcRotEnd;
                                              double fCalcRotWidth;
                                              double a2fCenter[2];
                                              bool   bIsIntersecting;
                                              
                                           
                                              if ((poRefln->bEllipsoidAvailable()) && (!poRefln->nPutGetEllipse(false,&a2x2fFindEllipseA[0][0],&a2fFindEllipseb[0],&fFindEllipsec,a2fCenter))) {
                                              
                                                  bEllipsoidInfoAvailable = true;
                                                  fFindRotMid = poRefln->fGetField(nFI2);
                                                  fFindRotWidth = poRefln->fGetField(nFI3);
                                                  fFindRotStart = fFindRotMid - 0.5*fFindRotWidth;
                                                  fFindRotEnd = fFindRotMid + 0.5*fFindRotWidth;
                                                  
                                                  fDSCOTH =  sqrt(fDStarSq  -  0.25f * fDStarSq*fDStarSq);
                                                  fCalcRotWidth = fLorentz*(fMosaicity[nComp][0] * fRADIANS_PER_DEGREE*fDSCOTH
                                                      + fMosaicity[nComp][1]*sqrt(fDStarSq-fDSCOTH*fDSCOTH))  + fMosaicity[nComp][2]*fRADIANS_PER_DEGREE;
                                                  

                                                  fCalcRotStart = fCalcRot - 0.5*fCalcRotWidth/Gs_dRADIANS_PER_DEGREE;
                                                  fCalcRotEnd = fCalcRot + 0.5*fCalcRotWidth/Gs_dRADIANS_PER_DEGREE;
                                                  bIsIntersecting = ((ABS(fFindRotEnd - fFindRotStart) + ABS(fCalcRotEnd - fCalcRotStart))>(max(fFindRotEnd,fCalcRotEnd) - min(fFindRotStart,fCalcRotStart)));
                                                  
                                                  if ((fEvalEllipse(fCalcPx0 - fObsPx0,fCalcPx1 - fObsPx1,&a2x2fFindEllipseA[0][0],&a2fFindEllipseb[0],&fFindEllipsec) <= 1.0) && bIsIntersecting) {
                                                      bEllipsoidIntersects = true;
                                                  }                                              
                                              }
                                          }
                                      }
                                  } else
                                      fRotErr = 0.0;
                              }                  

                              if ((fHKLErr<fMinTwinIDTolerance) && (fHKLErr<fBestHKLErr) && ((!bEllipsoidInfoAvailable) || (bEllipsoidIntersects)) ) {
                                  nBestTolerance = nComp;
                                  nBestTwinLaw = nTwinLaw;
                                  fBestRotErr = fRotErr;
                                  fBestHKLErr = fHKLErr;
                                  fBestMMErr = fMMErr;
                                  poRefln->vSetH(a3nHKL[0]);
                                  poRefln->vSetK(a3nHKL[1]);
                                  poRefln->vSetL(a3nHKL[2]);
                                  
                                  // If the packed HKL field is available, add the value.
                                  if (m_nFI_nPackedHKL >=0) 
                                      poRefln->vSetField(m_nFI_nPackedHKL,poRefln->nPackHKL());
                                  
                                  if ((m_nFI_fFloatH>=0) && (m_nFI_fFloatH>=0) && (m_nFI_fFloatH>=0)) {
                                      poRefln->vSetField(m_nFI_fFloatH,(float) a3fHKL[0]);
                                      poRefln->vSetField(m_nFI_fFloatK,(float) a3fHKL[1]);
                                      poRefln->vSetField(m_nFI_fFloatL,(float) a3fHKL[2]);
                                  }
                                  if (m_nFI_fHKLResid>=0) {
                                      poRefln->vSetField(m_nFI_fHKLResid,(float) fBestHKLErr);
                                  }
                              }
                          }                              
                      }
                      
                      // Fill in statistics.
                      if (nBestTolerance >= 0) {
                          afRMSMM[nBestTolerance][nBestTwinLaw] += fBestMMErr;
                          afRMSDEG[nBestTolerance][nBestTwinLaw] += fBestRotErr*fBestRotErr;
                          afRMSHKL[nBestTolerance][nBestTwinLaw] += fBestHKLErr;
                          if (poRefln->fGetSigmaI()>0.0) {
                              afIntensity[nBestTolerance][nBestTwinLaw] += poRefln->fGetIntensity();
                              afIoverSig[nBestTolerance][nBestTwinLaw] += poRefln->fGetIntensity()/max(1.0,poRefln->fGetSigmaI());
                              anIntensityContrib[nBestTolerance][nBestTwinLaw] ++;
                          }
                          anContrib[nBestTolerance][nBestTwinLaw]++;
                          nTotalContrib++;
                      }
                      

                      poRefln->vSetField(m_nFI_nTwinID,nBestTolerance+1 + 1000*nBestTwinLaw);
                      if (pfResidual)
                          pfResidual[i] = fBestHKLErr;
                  }
              }
         }     
     } else
         nStat = 1;

     if ((fRejectFraction>0.0) && (!nStat)) {
         // We need to reject a percentage of the data and recompute.
         qsort(&pfResidual[0],nRef,sizeof(float),float_cmp);
         fUserInput = fMinTwinIDTolerance = pfResidual[max(0,min(nRef-1,(int) ((nRef-1)*(1.0 - fRejectFraction))))];
         fRejectFraction = 0.0;
     } else
         break;
  } while (!nStat);

  // Update the error data.
  if (bComputeTwinID) {
      Ccrystal::nGetSetResidValues(oHeader,0,0,NULL,-1);
      for (nComp = 0; nComp < nComponents; nComp++) {
          for (nTwinLaw = 0; nTwinLaw <= anNumTwinLaws[nComp]; nTwinLaw++) {
              double afResidValues[g_nResidInfoRMSEntries];
              afResidValues[g_nResidInfoRMSMm] = sqrt(afRMSMM[nComp][nTwinLaw]/max(1,anContrib[nComp][nTwinLaw]));
              afResidValues[g_nResidInfoRMSHKL] = sqrt(afRMSHKL[nComp][nTwinLaw]/max(1,anContrib[nComp][nTwinLaw]));
              afResidValues[g_nResidInfoRMSRot] = sqrt(afRMSDEG[nComp][nTwinLaw]/max(1,anContrib[nComp][nTwinLaw]));
              afResidValues[g_nResidInfoNumIndexed] = anContrib[nComp][nTwinLaw];
              afResidValues[g_nResidInfoFraction] = anContrib[nComp][nTwinLaw]/((double) max(1,nTotalContrib));
              afResidValues[g_nResidInfoIntensity] = afIntensity[nComp][nTwinLaw]/max(1,anIntensityContrib[nComp][nTwinLaw]);
              afResidValues[g_nResidInfoIoverSig] = afIoverSig[nComp][nTwinLaw]/max(1,anIntensityContrib[nComp][nTwinLaw]);
              Ccrystal::nGetSetResidValues(oHeader,nComp + 1,nTwinLaw,&afResidValues[0],1);
          }
      }
  }
  


  // Delete all objects.
  delete m_poCrysGonio;
  delete m_poRotation;
  delete m_poSource;
  delete[] m_psDetectorNames;
  for (i=0;i<m_nNumDetectors;i++) {
      delete m_ppoDetector[i];
  }
  delete[] m_ppoDetector;
  if (pfResidual)
      delete[] pfResidual;

  return (nStat);
}


int
Creflnlist::nDeleteBatchID(Cstring& rsBatchIDTemplate)
{
  // Delete reflections in the list that match any in the rsBatchIDTemplate

  int nNum, nx, nStat;
  int *pnFlag;
  Cstring sTemp, sTemp2, sTemp3;
  Cstring sBatch;
  Crefln *poRefln;

  nNum = nGetNumReflns();
  pnFlag = new int [nNum];
  for (nx = 0; nx < nNum; nx++)
    pnFlag[nx] = 0;

  // Parse the rsBatchIDTemplate

  sTemp = rsBatchIDTemplate;
  //cout << "A. sTemp in nDeleteBatchID is: " << sTemp << endl << flush;
  while ("" != sTemp)
    {
      sTemp2 = sTemp.before(',');
      sTemp = sTemp.after(sTemp2);
      if (sTemp.GetAt(0) == ',') sTemp = sTemp.after(',');
      //cout << "B. sTemp in nDeleteBatchID is: " << sTemp << endl << flush;
      //cout << "sTemp2 in nDeleteBatchID is: " << sTemp2 << endl << flush;      
      if ("" != sTemp2)
        {
          if (sTemp2.contains('-'))
            {
              sTemp3 = sTemp2.after('-');  // No wildcards allowed!
              sTemp2 = sTemp2.before('-'); // No wildcards allowed!
              if ( sTemp2.contains('*') || sTemp2.contains('?') )
                {
                  cout << "ERROR, no wildcards allowed in a sBatch RANGE specification!\n"
                       << flush;
                  return (-1);
                }
              if ( sTemp3.contains('*') || sTemp3.contains('?') )
                {
                  cout << "ERROR, no wildcards allowed in a sBatch RANGE specification!\n"
                       << flush;
                  return (-1);
                }
              for (nx = 0; nx < nNum; nx++)
                {
                  poRefln = poGetRefln(nx);
                  sBatch = poRefln->sGetField(m_nFI_sBatch);
                  if ( (sTemp2 <= sBatch) && (sTemp3 >= sBatch) )
                    pnFlag[nx] = 1;
                }
            }
          else if ( sTemp2.contains('*') || sTemp2.contains('?') )
            {
              // Contains a * or ? wildcard

              for (nx = 0; nx < nNum; nx++)
                {
                  poRefln = poGetRefln(nx);
                  sBatch = poRefln->sGetField(m_nFI_sBatch);
                  if (sBatch.bRegularMatch(sTemp2))
                    pnFlag[nx] = 1;
                }
            }
          else 
            {
              // Does not contain -
              //cout << "B. sTemp2: " << sTemp2 << endl << flush;

              for (nx = 0; nx < nNum; nx++)
                {
                  poRefln = poGetRefln(nx);
                  sBatch = poRefln->sGetField(m_nFI_sBatch);
                  if (sTemp2 == sBatch)
                    pnFlag[nx] = 1;
                }
            }
        }
    }
    nStat = nDelete(1, pnFlag);
    //cout << "Batch reject: nNum, nStat: " << nNum << ", " << nStat << endl;
    return(0);
}

int Creflnlist::nDeleteRing(double* afIceRingReso,
                            int nRings, 
                            double fFull2ThetaWidthToReject, 
                            double fWavelength,
                            DTREK_WORD wCtrl)
{
  // Delete a group of reflections that occur in the same
  // resolution bins.  The idea is to remove reflns in ice or salt rings.
  // Current algorithm:
  // 1. Delete reflns in predefined rings

  // Old algorithm:
  // 1. Bin reflections into resolution bins.
  // 2. Find bins that have the most reflns.
  // 3. If bin contains 5% to 10% of all reflns, delete them.

  int nx;
  int nFI_fReso = -1;

  if( 0 <= m_nFI_fDetResolution && ! ( (wCtrl & DTREK_REFLNLIST_DL_RANK) && 
                                       (wCtrl & DTREK_REFLNLIST_DL_FROM_REFINE) 
                                     ) 
                                     )
    nFI_fReso = m_nFI_fDetResolution;
  else if (0 <= m_nFI_fResolution)
    nFI_fReso = m_nFI_fResolution;
  else
    return (0);  // No resolution field and we won't add one.

  if (fWavelength == 0.0)
      fWavelength = 1.0;

  float fReso;
  double fResoStart,fResoEnd;
  double f2ThetaMid;

  int nRing;
  int nCount = 0;
  int nRef;
  int *pnDelete = NULL;

  nRef = nGetNumReflns();
  Crefln *poRefln;
  pnDelete = new int [nRef];
  for (nx = 0; nx < nRef; nx++)   // Need to do this in case nRings = 0
    pnDelete[nx] = 0;

  for(nRing = 0; nRing < nRings; nRing++)
  {
      f2ThetaMid = 2.0 * asin(fWavelength/(2.0*afIceRingReso[nRing]));
      fResoStart = fWavelength / (2.0 * sin((f2ThetaMid - fFull2ThetaWidthToReject*Gs_dRADIANS_PER_DEGREE*0.5)/2.0));
      fResoEnd = fWavelength / (2.0 * sin((f2ThetaMid + fFull2ThetaWidthToReject*Gs_dRADIANS_PER_DEGREE*0.5)/2.0));
      
      for (nx = 0; nx < nRef; nx++)
      {
          if (nRing == 0)
            pnDelete[nx] = 0;

          poRefln      = poGetRefln(nx);
          fReso        = poRefln->fGetField(nFI_fReso);
          
          if( fReso <= fResoStart && fReso >= fResoEnd && !pnDelete[nx] )
          {
              pnDelete[nx] = 1;
              
              if( wCtrl & DTREK_REFLNLIST_DL_MARK_ONLY && m_nFI_nNonunfFlag > -1 )
                (*this)[nx].vSetField(m_nFI_nNonunfFlag, (int)enColorYellow);
   
              nCount++;
          }
      }
  }
  
  if( !(wCtrl & DTREK_REFLNLIST_DL_MARK_ONLY) )
    nDelete(1, pnDelete);
  
  delete [] pnDelete;
  pnDelete = NULL;

  if( !(wCtrl & DTREK_REFLNLIST_DL_MARK_ONLY) )
  {
    if( 0 < nCount )
        printf("\nINFO: %d reflns deleted out of %d that might be in ice rings.\n\n", nCount, nRef);
    else
        printf("\nINFO: No reflns deleted out of %d that might be in ice rings.\n\n", nRef);
  }

/********************************************************
Comment the following out for later use

  if (50 > m_nNumReflns)
    return (0);

  int nNumBins;
  int nBin;
  float afBins[101];
  float afBinCt[101];
  float fBinSize = 0.15;
  float fInt;

  if ("" != sGetEnv("DERING_BIN_SIZE"))
    nx = sscanf(sGetEnv("DERING_BIN_SIZE").string(), "%f", &fBinSize);
  if (1 != nx)
    fBinSize = 0.1;
  cout << "dering bin size is " << fBinSize << endl;
  for (nx = 0; nx < 101; nx++)
    {
      afBins[nx]  = 0.0;
      afBinCt[nx] = 0.0;
    }

  Crefln *poRefln;
  float fResoMax = 99999.0;
  float fResoMin = 0.0;
  float fReso;
  for (nx = 0; nx < m_nNumReflns; nx++)
    {
      poRefln  = *(m_ppoTheReflns+nx);
      fReso    = poRefln->fGetField(nFI_fReso);
      fResoMax = min(fResoMax, fReso);
      fResoMin = max(fResoMin, fReso);
      nBin     = (int)(fReso / fBinSize);
      if (nBin <   0) nBin =   0;
      if (nBin > 100) nBin = 100;
      afBins[nBin]  = afBins[nBin] + poRefln->fGetIntensity();
      afBinCt[nBin] = afBinCt[nBin] + 1.0;
    }

  // Now analyze the bins

  float fMinNum = 0.0;
  float fMaxNum = 0.0;

  fResoMin = 0.0;
  fResoMax = 0.0;
  for (nx = 0; nx < 101; nx++)
    {
      if (afBins[nx] > fMaxNum) 
        {
          fMaxNum = afBins[nx];
          fMinNum = (float)nx;
        }

      // Compute bin area normalizing factor as fReso and
      // multiply the number of reflns in a bin by the normalizing factor

      fReso    = (float)nx  * fBinSize;
      cout << nx << ", " << afBins[nx] << ", " << fReso << endl;
      if (0.0 < afBinCt[nx])
        afBins[nx] = afBins[nx] / afBinCt[nx];
      afBins[nx] = afBins[nx] * fReso;
      if (0.0 < afBins[nx])
        {
          fResoMin = fResoMin + afBins[nx];
          fResoMax = fResoMax + 1.0;
        }
    }
  fResoMin = fResoMin / fResoMax;
  cout << "Average density of reflns per bin-area: " << fResoMin << endl;

  cout << "Bin with max refln density: " << fMinNum
       << "  density refs: " << afBins[(int)fMinNum] << endl;

  fMinNum = fResoMin * 3.3; // Must be at least 3.3 times the average
  fMaxNum = (float)m_nNumReflns * 0.2; // 20%
  
  int *pnIndex;
  pnIndex = new int [m_nNumReflns];
  //for (nx = 0; nx < 101; nx++)
  float fLo, fHi;
  for (nx = 2; nx < 99; nx++)
    {
      fReso    = (float)nx  * fBinSize;      
      //if ( (nMinNum < afBins[nx]) && ((nMaxNum * fReso) > afBins[nx]) )
      //if (nMinNum < afBins[nx])
      
      fLo = ABS(afBins[nx-1]);
      if (fLo == 0.0) fLo = ABS(afBins[nx-2]);
      fHi = ABS(afBins[nx+1]);
      if (fHi == 0.0) fHi = ABS(afBins[nx+2]);
      fMaxNum = fLo + fHi;
      if (afBins[nx] > fMaxNum)
        {
          // Mark this reso bin for deletion
          
          fResoMin   = (float)nx * fBinSize;
          fResoMax   = fResoMin + fBinSize;
          cout << "Reso bin: " << fResoMax << " to " << fResoMin
               << " marked for deletion with density: " << afBins[nx]
               << endl;
          cout << "Lo, Hi, this_bin: " << fLo << ", " << fHi
               << ", " << afBins[nx];
          afBins[nx] = -afBins[nx];
        }
    }

  // Now mark reflns for deletion

  int nBins = 0;
  for (nx = 0; nx < m_nNumReflns; nx++)
    {
      poRefln     = *(m_ppoTheReflns+nx);
      fReso       = poRefln->fGetField(nFI_fReso);
      nBin        = (int)(fReso / fBinSize);
      pnIndex[nx] = 0;
      if (0.0 > afBins[nBin])
        {
          pnIndex[nx] = 1;
          nBins++;
        }
    }
  cout << "nDeleteRing() deleted " << nBins << " out of " 
       << m_nNumReflns << " reflns.\n\n";

  nDelete(1, pnIndex);
  delete [] pnIndex;
  pnIndex = NULL;
***********************************************************/
  return (0);
}





Cstring Creflnlist::ms_asTableHeaders[] = { " Resolution  \n    range",
                                            "     Int/sigmaI\n        range",
                                            "    Num\n reflns",
                                            "  Num\n rejs",
                                            "  I/sig\n    avg",
                                            "  Reso\n   max",
                                            " Backgnd\n     avg",
                                            " B/sig(B)\n      avg",
                                            " sharpness\n       avg",
                                            ""
                                           };

Cstring Creflnlist::ms_asTableFormat[] = { " %5.2lf -%5.2lf",
                                           " %c%4.0lf -  %c%4.0lf",
                                           " %6d",
                                           " %4d",
                                           " %6.1lf",
                                           " %5.2lf",
                                           " %7.2lf",
                                           " %8.1lf",
                                           " %9.2lf",
                                           ""
                                          };
   



struct tResoTableEntry
{
    double      m_fIoverSigAvg;
    double      m_fResoMax;
    int         m_nNumRefs;
    int         m_nNumRejs;
    double      m_fBackgroundAvg;
    double      m_fBackgroundOverSig;
    double      m_fSharpness;

    int operator ==(int nStat) { return 0; }
    int operator =(int nStat) { return 0; }
    
    void vClear() { memset(this,0,sizeof(*this)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////
int Creflnlist::nResoTable(double fResoMin, 
                           double fResoMax,
                           int nNumberOfResoBins,
                           int nNumberOfIOverSigmaBins,
                           int nTableTypes,
                           int nColsToPrint,
                           Cimage_header* poHeader,
                           DTREK_WORD wCtrl)
{
    if( nNumberOfResoBins < 1 || nNumberOfIOverSigmaBins < 1 )
        return 1; // safety

    itr<double>                 afValues;
    std::vector<Cstring>        asTableHeaders;
    std::vector<Cstring>        asTableFormats;
    int                         nx = 0;
    int                         nRef = 0;
    
    bool                        bForRanking = false;
    Cstring     strCallingModuleName = sDtrekGetModuleName();
    strCallingModuleName.upcase();

    if( strCallingModuleName == "DTFIND" )
        bForRanking = true;
    else if( strCallingModuleName == "DTREFINE" ) 
        bForRanking = (wCtrl & DTREK_REFLNLIST_RT_RANK) > 0U ? true : false;
    else 
        return 1; // for now unsupported calling module

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Decide on where to get the resolution from...
    int         nFI_fReso = 0;
    if( strCallingModuleName == "DTFIND" || (strCallingModuleName == "DTREFINE" && !bForRanking) )
        nFI_fReso = m_nFI_fDetResolution >= 0 ? m_nFI_fDetResolution : m_nFI_fResolution;
    else
        nFI_fReso = m_nFI_fResolution;
    ///////////////////////////////////////////////////////////////////////////////////////////////

    double          f0 = 0.0, f1 = 0.0;
    Cstring         sTemp("");
    
    if( fResoMin < 0.0 || fResoMax < 0.0 )
    {
        // Determine the resolution from the reflections.
        if( nFI_fReso >= 0 )
            vGetResoLimits(fResoMin, fResoMax, nFI_fReso);
        else
            return 1;
    }
    
    if( fResoMin < 0.0 || fResoMax < 0.0 )
        return 1;

    if( fResoMin < fResoMax )
      std::swap(fResoMin,fResoMax);
    
    // Build the I/sig bins and/or resolution bins.


    itr<double> afResoStart;
    itr<double> afResoEnd;
    itr<double> afIoverSigStart;
    itr<double> afIoverSigEnd;
    tResoTableEntry oAverageStats;
    itr<tResoTableEntry> aoResoStats;
    itr<tResoTableEntry> aoIoverSigStats;
    
    const int   nIoverSigBins = 10;
    bool bIsDtRefine = false;


    // Determine resolution bins.
    f0 = 1.0/pow(fResoMin,3.0);   
    f1 = (1.0/pow(fResoMax,3.0) - 1.0/pow(fResoMin,3.0))/nNumberOfResoBins; 
    for (nx = 0; nx < nNumberOfResoBins; nx++) {
        afResoStart + min(99.99,(1.0/ max(1e-10,pow(f0,0.3333333))));
        afResoEnd + (1.0/ max(1e-10,pow(f0 + f1,0.3333333)));
        f0 += f1;
    }
      

    // Determine I/sig bins.
    -afValues;
    float       fIntensity = 0.0;
    float       fSigma = 0.0;
    for (nRef = 0; nRef <  m_nNumReflns; nRef++)
    {
        fIntensity = (*this)[nRef].fGetIntensity();
        fSigma = (*this)[nRef].fGetSigmaI();
        
        if( bForRanking && 0.0 >= fIntensity )  //RB: get rid of "bad" reflections ( -999.000 == fIntensity && -999.000 == fSigma )
            continue;
        
        afValues + fIntensity / max(1.0, fSigma);
    }
    
    if( afValues.size() > 0 ) // safety check 
        qsort(&afValues[0], afValues.size(), sizeof(double), double_cmp);
    
    f0 = max(2,(int) ceil(afValues[(int) (afValues.size()/2.0)]/(nIoverSigBins - 1)));
    
    -afIoverSigStart;
    -afIoverSigEnd;
    
    for (nx = nIoverSigBins - 1; nx >= 0; nx--) {
        afIoverSigStart + (nx*f0);
        afIoverSigEnd + ((nx+1)*f0);
    }

    // Initialize table entries.
    aoResoStats.setsize(afResoStart.size());
    aoIoverSigStats.setsize(afIoverSigStart.size());
    for (nx = 0; nx < aoResoStats.size(); nx++)
        aoResoStats[nx].vClear();
    for (nx = 0; nx < aoIoverSigStats.size(); nx++)
        aoIoverSigStats[nx].vClear();
    oAverageStats.vClear();

    // The main computing loop.

    double              fIoverSig = 0.0;
    double              fReso = 0.0;
    int                 nIoverSig = 0;
    int                 nReso = 0;
    int                 nTable = 0;
    tResoTableEntry*    poStats = NULL;
    
    for (nRef = 0; nRef < m_nNumReflns; nRef++) 
    {
        /////////////////////////////////////////////////////////////////////////
        fIntensity = (*this)[nRef].fGetIntensity();
        fSigma = (*this)[nRef].fGetSigmaI();
        
        if( bForRanking && 0.0 >= fIntensity )  //RB: get rid of "bad" reflections ( -999.000 == fIntensity && -999.000 == fSigma )
           continue;
        /////////////////////////////////////////////////////////////////////////

        fReso = (*this)[nRef].fGetField(nFI_fReso);
        fIoverSig = fIntensity / max(1.0, fSigma);

        for (nReso = 0; (afResoEnd[nReso]>fReso) && (nReso + 1 < afResoEnd.size()); nReso++);
        for (nIoverSig = 0; (afIoverSigEnd[nIoverSig] > fIoverSig) && (nIoverSig + 1 < afIoverSigEnd.size()); nIoverSig++);

        
        for(nTable = 0; nTable < 3; nTable ++)
        {
            if (nTable == 0) 
                poStats = &aoResoStats[nReso];
            else if (nTable == 1)
                poStats = &aoIoverSigStats[nIoverSig];
            else
                poStats = &oAverageStats;
            
            if( (*this)[nRef].m_nSelect == m_pcSelect[m_nTotalFieldsPlus1-1] )
            {
                if (!poStats->m_nNumRefs) 
                    poStats->m_fResoMax = fReso;
                else
                    poStats->m_fResoMax = min(fReso,poStats->m_fResoMax);
                poStats->m_fIoverSigAvg += fIoverSig;
                
                if ((m_nFI_fBackgroundSigma>=0) && (m_nFI_fBackground>=0))
                    poStats->m_fBackgroundOverSig += (*this)[nRef].fGetField(m_nFI_fBackground)/max(0.01,(*this)[nRef].fGetField(m_nFI_fBackgroundSigma));
                
                if (m_nFI_fBackground>=0)
                    poStats->m_fBackgroundAvg += (*this)[nRef].fGetField(m_nFI_fBackground);
                
                if (m_nFI_fObsSharpness)
                    poStats->m_fSharpness += (*this)[nRef].fGetField(m_nFI_fObsSharpness);

                poStats->m_nNumRefs++;
                
            }
            else 
            {
                poStats->m_nNumRejs++;
                bIsDtRefine = true;
            }
        }

    }

    // Now, build the tables for printing.
    int nRow;
    int nRows;
    int nTableType;

    tResoTableEntry* poRow;

    for (nTableType = 1; nTableType <=2; nTableType = nTableType << 1) 
    {
        if (nTableTypes & nTableType)
        {
            if (nTableType & eTableVsResolution)
            {
                nRows = aoResoStats.size();
                poRow = & aoResoStats[0];
                if (nColsToPrint == -1)
                {
                    nColsToPrint = 
                        eTableResoStartEnd | eTableNumReflns | eTableNumRejects | eTableIoverSig | eTableAvgBackground | 
                        eTableAvgBackgroundOverAvgBackgroundSig | eTableSharpness;
                }
            } 
            else if (nTableType & eTableVsIoverSig)
            {
                nRows = aoIoverSigStats.size();
                poRow = & aoIoverSigStats[0];
                if (nColsToPrint == -1)
                {
                    nColsToPrint = 
                        eTableIoverSigRange | eTableNumReflns | eTableNumRejects | eTableMaxReso | eTableAvgBackground | 
                        eTableAvgBackgroundOverAvgBackgroundSig | eTableSharpness;
                }
            } 
            else
                continue;
            
            asTableHeaders.clear();
            asTableFormats.clear();
            for (nx = 0; ms_asTableFormat[nx].length();nx++)
            {
                if ((nColsToPrint >> nx) & 1)
                {
                    asTableHeaders.push_back((Cstring) ms_asTableHeaders[nx]);
                    asTableFormats.push_back((Cstring) ms_asTableFormat[nx]);
                }
            }

            -afValues;
            double      fAverageSharpness = 0.0;
            
            for (nRow = 0; nRow <= nRows; nRow++,poRow++)
            {
                if (nRow == nRows)
                    poRow = &oAverageStats;
                
                if (nColsToPrint & eTableResoStartEnd)
                {
                    if (nRow == nRows)
                        afValues + g_fTableMax + g_fTableMin;
                    else
                        afValues + afResoStart[nRow] + afResoEnd[nRow];
                }
                
                if (nColsToPrint & eTableIoverSigRange)
                {
                    if (nRow == nRows)
                        afValues + ((double) '<') + g_fTableMin + ((double) '>') + g_fTableMax;
                    else
                        afValues + ((double) ((nRow == nRows-1)?'<':' ')) + afIoverSigStart[nRow] + ((double) ((nRow == 0)?'>':' ')) + afIoverSigEnd[nRow];
                }
                
                if (nColsToPrint & eTableNumReflns)
                {
                    afValues + poRow->m_nNumRefs;
                }
                
                if (nColsToPrint & eTableNumRejects)
                {
                    afValues + poRow->m_nNumRejs;
                }
                
                if (nColsToPrint & eTableIoverSig) 
                {
                    afValues + poRow->m_fIoverSigAvg/max(1,poRow->m_nNumRefs);
                }
                
                if (nColsToPrint & eTableMaxReso)
                {
                    afValues + poRow->m_fResoMax;
                }
                
                if (nColsToPrint & eTableAvgBackground)
                {
                    afValues + poRow->m_fBackgroundAvg/max(1,poRow->m_nNumRefs);
                }
                
                if (nColsToPrint & eTableAvgBackgroundOverAvgBackgroundSig)
                {
                    afValues + poRow->m_fBackgroundOverSig/max(1,poRow->m_nNumRefs);
                }
                
                if (nColsToPrint & eTableSharpness)
                {
                    afValues + poRow->m_fSharpness/max(1,poRow->m_nNumRefs);
                    fAverageSharpness += afValues.last();
                }
            }

            sTemp = (Cstring) ((nTableType & eTableVsResolution)?"Resolution statistics":"Intensity statistics");
            
            nPrintTable(sTemp, asTableHeaders, asTableFormats, afValues, true);
            
            nColsToPrint = -1;
        }
    }
    
    return 0;
}

void Creflnlist::vObliqueIncidenceCorrection(const int nApply, Cimage_header *poHeader)
{
  // Most of this routine was taken from Cintegrate.cpp

  // Apply an oblique incidence correction to the fIntensity and fSigmaI
  // fields.  The reflnlist should be a dtprofit.ref or a dtintegrate.ref
  // that has the fCalc_oblique field in it
  // If nApply = 0, then APPLY the correction.  This is the default
  // If nApply <> 0, then UNAPPLY the correction.

  if (0 > m_nFI_fOblique)
    {
      cout << "ERROR - no " << ms_sfOblique << " field in the reflnlist.\n";
      return;
    }
  if (NULL == poHeader)
    {
      cout << "ERROR - no header info available.\n";
      return;
    }
  if (!poHeader->bIsAvailable())
    {
      cout << "ERROR - no header info available.\n";
      return;
    }

  // Apply any oblique incidence correction
  // See for example: Zaleski, Wu & Coppens (1998)
  //                     J. Appl. Cryst. 31, 302-304.
  // Icorr = Iobs * Kobl;
  // Tobl  = (1.0 - exp(Foblnorm)) / (1.0 - exp(Foblnorm/cosA))
  // The refln has 1/cosA in the fCalc_oblique field.
  // m_nNumObliqueCorr may equal 0.
  // Remember (1.0 - exp(Foblnorm)) = m_pfObliqueCorrNorm
  
  // !!!Note: The above is a TRANSMISSION effect.
  //          An ABSORPTION effect would be
  // Aobl = exp(Foblnorm) / exp(Foblnorm/cosA)
  //      = -exp(Foblnorm) / -exp(Foblnorm/cosA)
  
  int nx, i;
  int nObl;
  int nStat;
  int nNumObliqueCorr;
  Crefln *poRefln;
  float fRefOblique;
  float *pfObliqueCorr;
  float *pfObliqueCorrNorm;
  float *pfObliqueCorrType;
  float *pfTemp;

  Cstring sTemp = "";
  nNumObliqueCorr = 0;
  nStat = poHeader->nGetValue(D_K_IntegrateOblique, &nNumObliqueCorr);
  if (0 >= nNumObliqueCorr)
    {
      cout << "INFO: no oblique incidence info in input header.\n";
      return;
    }

  pfObliqueCorr     = new float [nNumObliqueCorr+1];
  pfObliqueCorrNorm = new float [nNumObliqueCorr+1];
  pfObliqueCorrType = new float [nNumObliqueCorr+1];

  pfTemp = new float [2*nNumObliqueCorr + 1];
  nStat = poHeader->nGetValue(D_K_IntegrateOblique, 2*nNumObliqueCorr+1,
                              pfTemp);
  int nOff = 1;
  for (i = 0; i < nNumObliqueCorr; i++)
    {
      // Shift the values, since the first is really the number of factors

      // 0 means absorp, 1 means Trans

      pfObliqueCorrType[i] = pfTemp[nOff++];
      pfObliqueCorr[i]     = pfTemp[nOff++];

      if (   (0.0 != pfObliqueCorrType[i])
          && (1.0 != pfObliqueCorrType[i]))
        {
          cout << "WARNING, oblique incidence correction type is not 0 or 1!" << endl;
        }

      // Calculate 1 - exp(factor), when angle is normal to detector plane.
      // for Transmission types
      // Calculate -exp(factor), when angle is normal to detector plane.
      // for Absorption types

      pfObliqueCorrNorm[i] = pfObliqueCorrType[i]
                                   - (float)exp(pfObliqueCorr[i]);
      if(0.0f == pfObliqueCorrNorm[i])
        {
          cout << "ERROR: Oblique incidence correction factor SHOULD NOT be 0!\n";
          if (1 == nNumObliqueCorr)
            {
              nNumObliqueCorr = 0;
              cout << "       NOT USING any OBLIQUE INCIDENCE CORRECTION!\n\n";
            }
          else
            {
              cout << "       Reset to normal to detector factor to 1.0!\n"
                   << "       ALL RESULTS WILL BE BOGUS!\n\n";
              pfObliqueCorrNorm[i] = 1.0f;
            }
        }
    }
  delete [] pfTemp;
  pfTemp = NULL;
  cout << "\nINFO: Num oblique incidence factors:" << nNumObliqueCorr;
  for (i = 0; i < nNumObliqueCorr; i++)
    {
      cout << "\n    Oblique incidence factor " << i+1 << ":   "
           << pfObliqueCorr[i] <<  "   Norm: "
           << pfObliqueCorrNorm[i]
           << "   Type: "
           << pfObliqueCorrType[i];
    }
  if (0 == nApply)
    cout << "\nThe above factors will be APPLIED.\n";
  else
    cout << "\nThe above factors will be UNAPPLIED.\n";

  for (nx = 0; nx < m_nNumReflns; nx++)
    {
      poRefln     = *(m_ppoTheReflns+nx);
      fRefOblique = 1.0f;
      for (nObl = 0; nObl < nNumObliqueCorr; nObl++)
        {
          fRefOblique *= (pfObliqueCorrNorm[nObl]
                        / (pfObliqueCorrType[nObl]
                   - (float)exp(pfObliqueCorr[nObl]
        * poRefln->fGetField(m_nFI_fOblique))));
        }
      if (0 == nApply)
        {
          poRefln->vSetIntensity(poRefln->fGetIntensity() *fRefOblique);
          poRefln->vSetSigmaI(poRefln->fGetSigmaI() * fRefOblique);
        }
      else
        {
          // Unapply the correction

          poRefln->vSetIntensity(poRefln->fGetIntensity() / fRefOblique);
          poRefln->vSetSigmaI(poRefln->fGetSigmaI() / fRefOblique);
        }
    }
  delete [] pfObliqueCorr; pfObliqueCorr = NULL;
  delete [] pfObliqueCorrNorm; pfObliqueCorrNorm = NULL;
  delete [] pfObliqueCorrType; pfObliqueCorrType = NULL;
  return;
}

void Creflnlist::vAnalyzeBkgVsResolution()
{
// RB Disabled, as it is still under development...
#if 0
    int     nReflnCount = nGetNumReflns();

    int*    pnIndex = (int*)malloc(sizeof(int)*nReflnCount);

    //    vSort(eReflnField_float_type, m_nFI_fResolution, pnIndex);
    vSort(eReflnField_float_type, m_nFI_fBackground, pnIndex);

    //   nWrite(sDtrekGetPrefix() + "dtrefpredbkg.ref");
    
    static      int     nFileCount = 0;
    
    Cstring     strFileCount(++nFileCount);
    Cstring     strFileName = "BKG_" + strFileCount;
    strFileName += ".ref";
    
    FILE*   pFile = fopen(strFileName, "w");
    char    szOutLine[256];
    
    //float   fRes = 0.0;
    float   fBkg = 0.0;  
    float   fSigBkg  = 0.0;
    
    int     nH = 0;
    int     nK = 0;
    int     nL = 0;

    for(int ii=0; ii < nReflnCount; ii++)
    {
        Crefln*     pRefln = poGetRefln(pnIndex[ii]); 

        //fRes = pRefln->fGetField(m_nFI_fResolution);
        fBkg = pRefln->fGetField(m_nFI_fBackground);
        fSigBkg = pRefln->fGetField(m_nFI_fBackgroundSigma);
        
        if( fBkg < 0.0 || fSigBkg < 0.0 )
            continue;

        sprintf(szOutLine, "%d %d %d %.5f %.5f %.5f %.5f %.5f", 
                pRefln->nGetH(),
                pRefln->nGetK(),
                pRefln->nGetL(),
                pRefln->fGetField(m_nFI_fResolution),
                fBkg, 
                fSigBkg, 
                log((double)fBkg)/log((double)fSigBkg), 
                (double)fSigBkg/sqrt((double)fBkg));
        
        if( ii != nReflnCount-1 )
            strcat(szOutLine, "\n");

        fprintf(pFile, "%s", szOutLine);
    }

    fclose(pFile);
    free (pnIndex);
#endif
}
void Creflnlist::vSetExternalSortIndex(int* pnSortIndex)
{
    if (NULL != m_pnSortIndex)
    {
        delete [] m_pnSortIndex;
        m_pnSortIndex = NULL;
    }
    
    m_pnSortIndex = pnSortIndex;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Loop over the reflection list, do statistics on ranking parameters and output that info into the header
void Creflnlist::vRank(Cimage_header* poHeader,
                       double fResoMin,      
                       double fResoMax,      
                       int nNumberOfResoBins,
                       int nMinReflnsForStats,
                       bool bEvalExpTime)
{
    if( !poHeader )
        return; // safety
    
    // Figure out the calling module name and the field number for the resolution.
    Cstring     strCallingModuleName = sDtrekGetModuleName();
    strCallingModuleName.upcase();
    
    // Decide on where to get the resolution from...
    int         nFI_fReso = 0;
    if( strCallingModuleName == "DTFIND" )
        nFI_fReso = m_nFI_fDetResolution >= 0 ? m_nFI_fDetResolution : m_nFI_fResolution;
    else
        nFI_fReso = m_nFI_fResolution;
    
    strCallingModuleName += "_RANK";

    // Do we have valid resolution limits?
    if( fResoMin < 0.0 || fResoMax < 0.0 )
    {
        // Determine the resolution from the reflections.
        if( nFI_fReso >= 0 )
            vGetResoLimits(fResoMin, fResoMax, nFI_fReso);
        else
        {
            printf("ERROR: Failed to rank reflection list. Resolution field is invalid.\n");
            return;
        }
    }
    
    if( fResoMin < 0.0 || fResoMax < 0.0 )
    {
        printf("ERROR: Failed to rank reflection list. Invalid resolution limits.\n");
        return;
    }

    if( fResoMin < fResoMax )
      std::swap(fResoMin,fResoMax);

    ////////////////////////////
    // Generate resolution bins
    CResoBins<CResoRankingBin>  oRankingResoBins(fResoMin, fResoMax, nNumberOfResoBins);
    
    if( !oRankingResoBins.bIsAvailable() )
    {
        printf("ERROR: Failed to rank reflection list. Cannot set up resolution bins.\n");

        return;
    }
    
    CResoRankingBin*            pResoBin = NULL;

    for(int ii=0; ii < nNumberOfResoBins; ii++)
    {
        pResoBin = oRankingResoBins.pGetBin(ii);
        
        pResoBin->vSetMinNumberOfReflnsForStats(nMinReflnsForStats);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Generate rotation bins. This is needed for ranking anisotropy analysis.
    // We assume that we need at least 10 degrees of rotation to detect anisotropy.
    const   double      c_dMinRot  = -360.0;   // degrees, the lowest possible rotation
    const   double      c_dMaxRot  =  360.0;   // degrees, the largest possible rotation
    const   double      c_dRotStep =   10.0;   // degrees, chosen granularity for the bin statistics 
    std::vector<int>        anRotBinFlag;      // flag=0 means: there are no reflns in that bin; flag=1 means there are some.
    double      dCur = c_dMinRot;
    while(dCur <= c_dMaxRot)
    {
        anRotBinFlag.push_back(0);  // Initialize each rotation bin as "inactive", i.e. not having any reflns

        dCur += c_dRotStep;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    double    dIntensity               = 0.0;
    double    dSigma                   = 0.0;
    double    dPeakArea                = 0.0;
    double    dPeakBorderCount         = 0.0;
    
    double    dMidRotation             = 0.0;
    

    double    dMajorAxis               = 0.0;
    double    dMinorAxis               = 0.0;

    double    dEllipsArea              = 0.0;

    double    dProfileAsymmetry        = 0.0;
    double    dMeanProfileAsymmetry    = 0.0;
    
    double    dMeanIntDensity          = 0.0;
    
    double    dMeanShapeIrregularity   = 0.0;
    double    dShapeIrregularity       = 0.0;
    double    dEffectivePeakAreaRadius = 0.0;

    double    dMeanShapeAnisotropy     = 0.0;

    int       nUsedReflns              = 0;

    double    dReso  				   = 0.0;

    int                         nRotBinIndex = 0;

    for(int nRef = 0; nRef < m_nNumReflns; nRef++) 
    {
        dIntensity = (*this)[nRef].fGetIntensity();
        dSigma     = (*this)[nRef].fGetSigmaI();

        if( 0.0 >= dIntensity || 0.0 >= dSigma )  //RB: for ranking purposes get rid of "bad" reflections
           continue;

        if( (*this)[nRef].m_nSelect != m_pcSelect[m_nTotalFieldsPlus1-1] )
            continue;            // RB: not sure why Thad needed this
        
        if( -1 < m_nFI_fObsRotMid && -1 < nFI_fReso )
        {
            dMidRotation = (double)((*this)[nRef].fGetField(m_nFI_fObsRotMid));

            if( dMidRotation >= c_dMinRot && dMidRotation <= c_dMaxRot )
            {
                // Figure out the rotation bin index and mark that rotation bin as active
                nRotBinIndex = (int)floor((dMidRotation - c_dMinRot) / c_dRotStep);
                anRotBinFlag[nRotBinIndex] = 1;
            
                ////////////////////////////////////////////
                // Get resolution
                dReso = (*this)[nRef].fGetField(nFI_fReso);
            
                ////////////////////////////////////////////
                // Add IOverSigmaValue to resolution bin
                pResoBin = oRankingResoBins.pGetBin(1.0 / pow((double)dReso, 3.));
            
                if( pResoBin )
                {
                    pResoBin->vAddIOverSigmaValue(dMidRotation, dIntensity / dSigma); // dSigma is already tested for zero above
                    pResoBin->vAddToIntensityStats(dIntensity);
                    pResoBin->vAddToIOverSigmaStats(dIntensity/dSigma);
                } 
           }
        }

        if( -1 < m_nFI_fObsSharpness )
        {
            dProfileAsymmetry = (*this)[nRef].fGetField(m_nFI_fObsSharpness);
            dMeanProfileAsymmetry += dProfileAsymmetry;
        }
        
        if( -1 < m_nFI_fObsPeakAreaCount )
            dPeakArea = (*this)[nRef].fGetField(m_nFI_fObsPeakAreaCount);
                 
        if( -1 < m_nFI_fObsPeakBorderCount )
            dPeakBorderCount = (*this)[nRef].fGetField(m_nFI_fObsPeakBorderCount);

        if( dPeakArea > 0.0 )
        {
            dEffectivePeakAreaRadius = sqrt(dPeakArea / Gs_dPI);
            
            // Shape amorphousness is characterized by the circumference/area ratio.
            // The "ideal" spot is a circle of radius r, so the ideal ratio = 2 * PI * r / PI * r^2 = 2 / r.
            // So we will take the ratio of the *real* ratio to the *ideal* one and compare that to 1.0 
            // 10 is just an arbitrary scale factor. The *real* ratio below is dPeakBorderCount/dPeakArea.
            
            dShapeIrregularity = 10.0 * ( dPeakBorderCount * dEffectivePeakAreaRadius / 2.0 / dPeakArea  - 1.0 );  
            
            if( dShapeIrregularity < 0.0 )  // it shouldn't really!
                dShapeIrregularity = 0.0;

            dMeanShapeIrregularity += dShapeIrregularity;
        }

        // Do ellips numbers
        dMajorAxis             = 0.0; 
        dMinorAxis             = 0.0;         
        dEllipsArea            = 0.0;         

        //double dEffEllipsCircumf       = 0.0;
        //double dEffEllipsArea          = 0.0;
        //double dEffEllipsCircumfToArea = 0.0;
        //double dEffEllipsMajor         = 0.0;
        //double dEffEllipsMinor         = 0.0;

        if( 0 == (*this)[nRef].nConvertToAxialEllipse(&dMajorAxis, &dMinorAxis) )
        {
            if( dMajorAxis < dMinorAxis )
	      std::swap(dMajorAxis, dMinorAxis);
        
            if( 0.0 != dMinorAxis )
            {
                dMeanShapeAnisotropy += (dMajorAxis / dMinorAxis - 1.0); // -1.0, because a circle has the ratio equal to 1.0
            }

            dEllipsArea = dMajorAxis * dMinorAxis / 4.0 * fPI;
        }
        
        if( dEllipsArea > 0.0 )
        {
            dMeanIntDensity += dIntensity / dEllipsArea;
         }

        //        if( dPeakArea > 0.0 )
        //        {
        //            
        //            dEffEllipsMinor =  sqrt(dPeakArea / Gs_dPI / (dMajorAxis / dMinorAxis));
        //            dEffEllipsMajor =  dEffEllipsMinor * (dMajorAxis / dMinorAxis);
        //            dEffEllipsArea = Gs_dPI * dEffEllipsMajor * dEffEllipsMinor;
        //            dEffEllipsCircumf = Gs_dPI * ( 1.5 * (dEffEllipsMajor + dEffEllipsMinor) 
        //                                - sqrt(dEffEllipsMajor*dEffEllipsMinor) );
        //            dEffEllipsCircumfToArea = dEffEllipsCircumf / dEffEllipsArea;
        //
        //            dShapeIrregularity = 10 * ( dPeakBorderCount / dPeakArea / dEffEllipsCircumfToArea  - 1.0 );  
        //
        //            if( dShapeIrregularity < 0.0 )  // it shouldn't really!
        //                dShapeIrregularity = 0.0;
        //
        //            dMeanShapeIrregularity += dShapeIrregularity;
        //        }
        
        nUsedReflns++;
    }
    
    if( 0 == nUsedReflns )
        return;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    Cstring     strTemp("");
    Cstring     strParamType("");
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    dMeanProfileAsymmetry       /= nUsedReflns;
    //oRotBinsProfileAsymmetry.bGetStatistics(dRotationBinMean, dRotationBinDev);
    
    strTemp = dMeanProfileAsymmetry;
    //strTemp += ' ';
    //strTemp += dRotationBinDev;

    strParamType = "PEAK_PROFILE_ASYMMETRY";
    poHeader->nReplaceValueDictionary(strCallingModuleName, strParamType, strTemp);
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    dMeanIntDensity             /= nUsedReflns;
    //oRotBinsIntDensity.bGetStatistics(dRotationBinMean, dRotationBinDev);
    
    strTemp = dMeanIntDensity;
    //strTemp += ' ';
    //strTemp += dRotationBinDev;
    
    strParamType = "PEAK_INTENSITY_DENSITY";
    poHeader->nReplaceValueDictionary(strCallingModuleName, strParamType, strTemp);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    dMeanShapeIrregularity      /= nUsedReflns;
    //oRotBinsShapeIrregularity.bGetStatistics(dRotationBinMean, dRotationBinDev);
    
    strTemp = dMeanShapeIrregularity;
    //strTemp += ' ';
    //strTemp += dRotationBinDev;
    
    strParamType = "PEAK_SHAPE_IRREGULARITY";
    poHeader->nReplaceValueDictionary(strCallingModuleName, strParamType, strTemp);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    dMeanShapeAnisotropy        /= nUsedReflns;
    //oRotBinsShapeAnisotropy.bGetStatistics(dRotationBinMean, dRotationBinDev);
    
    strTemp = dMeanShapeAnisotropy;
    //strTemp += ' ';
    //strTemp += dRotationBinDev;
    
    strParamType = "PEAK_SHAPE_ANISOTROPY";
    poHeader->nReplaceValueDictionary(strCallingModuleName, strParamType, strTemp);
    /////////////////////////////////////////////////////////////////////////////////

    Cstring     strResoShellKeyWord = "RESO_SHELL_*";
    Cstring     strHeaderValue("");
    char        szTemp[256];

    double      dIOverSigmaOverallMean = -1.0;
    double      dIOverSigmaRotationDev = -1.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Based on the vector of rotation bins flags, generate a vector of active bins
    double          dBinStart = -1.0;
    double          dBinEnd   = -1.0;
    std::vector<CSegment>     vecActiveRotBins;
    
    for(int iRotBin=0; iRotBin < anRotBinFlag.size(); iRotBin++)
    {
        if( 0 != anRotBinFlag[iRotBin] )
        {
            dBinStart = c_dMinRot + iRotBin * c_dRotStep;
            dBinEnd   = dBinStart + c_dRotStep;

            vecActiveRotBins.push_back( CSegment(dBinStart, dBinEnd) );
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Output some general information about resolution bins into the header
    strHeaderValue = nNumberOfResoBins;
    strHeaderValue += ' ';
    strHeaderValue += Cstring(nMinReflnsForStats);
    poHeader->nReplaceValueDictionary(strCallingModuleName, "RESO_SHELLS", strHeaderValue);
    strHeaderValue = "";
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Loop over all reso bins, and for every bin calculate a mean I/Sigma value and its deviation over rotation ranges 
    // and output that information into the header.
    for(int iResoBin=0; iResoBin < nNumberOfResoBins; iResoBin++)
    {
        pResoBin = oRankingResoBins.pGetBin(iResoBin);
        
        pResoBin->vCalcIOverSigma(vecActiveRotBins);

        dIOverSigmaOverallMean = pResoBin->fGetIOverSigmaMean();
        dIOverSigmaRotationDev = pResoBin->fGetIOverSigmaMeanStandardDeviationWRTRotation();

        sprintf(szTemp, "%.2f %.2f", pResoBin->dGetLowReso(), pResoBin->dGetHighReso());
        strHeaderValue = szTemp;

        sprintf(szTemp, " %d", pResoBin->nGetReflnCount());
        strHeaderValue += szTemp;

        sprintf(szTemp, " %.1f %.1f", dIOverSigmaOverallMean, dIOverSigmaRotationDev);
        strHeaderValue += szTemp;

        poHeader->nReplaceValueDictionary(strCallingModuleName, strResoShellKeyWord, strHeaderValue, iResoBin);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Do a linear fit of log(I/Sigma) vs. (sin(Theta)/lambda)^2
// Using the fit results, find I/Sigma value for the target resolution. 
// Knowing the image exposure time and how I/Sigma behaves as a function of exposure time,
// figure out what exposure time is needed to bring I/Sigma to the target level.
bool Creflnlist::bEvaluateExpTime(Cimage_header* poHeader, 
                                  double fTargetReso, 
                                  double fTargetIOverSigma, 
                                  double fIOverSigmaTimePower)
{
    if( !poHeader || fTargetReso < 0.0 )
        return false; // safety
    
    // Figure out the calling module name and the field number for the resolution.
    Cstring     strCallingModuleName = sDtrekGetModuleName();
    strCallingModuleName.upcase();
    
    // Decide on where to get the resolution from...
    int         nFI_fReso = 0;
    if( strCallingModuleName == "DTFIND" )
        nFI_fReso = m_nFI_fDetResolution >= 0 ? m_nFI_fDetResolution : m_nFI_fResolution;
    else
        nFI_fReso = m_nFI_fResolution;
    if( 0 > nFI_fReso )
        return false;
    
    strCallingModuleName += "_RANK";
    
    std::vector<float>      afLogIOverSigmaValues;
    std::vector<float>      afWeights;
    std::vector<float>      afSsquared;
    
    // Since numerical recipes require an array to start from index 1, we need to fill index 0 with a value.
    afLogIOverSigmaValues.push_back(-1.0f);
    afWeights.push_back(-1.0f);
    afSsquared.push_back(-1.0f);
    
    double      fReso = -1.0;
    
    double      fIntensity = -1.0;
    double      fSigma = -1.0;
    double      fLogIOverSigmaValue = -1.0;

    for(int nRef = 0; nRef < m_nNumReflns; nRef++) 
    {
        if( (*this)[nRef].m_nSelect != m_pcSelect[m_nTotalFieldsPlus1-1] )
            continue;            // RB: not sure why Thad needed this
        
        fIntensity = (*this)[nRef].fGetIntensity();
        fSigma     = (*this)[nRef].fGetSigmaI();

        if( 0.0 >= fIntensity || 0.0 >= fSigma )  //get rid of "bad" reflections
           continue;

        fReso = (*this)[nRef].fGetField(nFI_fReso);
        if( fReso <= 0.0 )// just a safety check
            continue;
        else if( fReso > cfMinimumResolutionForLinearApproximation )
            continue;  // lower bins do not exhibit linear dependency of log(I/Sigma(I)) vs resolution

        fLogIOverSigmaValue = fIntensity / fSigma;
        
        afLogIOverSigmaValues.push_back((float)log10(fLogIOverSigmaValue));
        afWeights.push_back(1.0f);

        afSsquared.push_back(0.25f/(float)fReso/(float)fReso);
    }
    
    if( afLogIOverSigmaValues.size() - 1 < cnMinimumNumberOfPointsForFitting )
        return false;  // no sense fitting a straight line if we have too few points
    
    float       fA = 0.0f;
    float       fB = 0.0f;
    float       fSigA = 0.0f;
    float       fSigB = 0.0f;
    float       fChiSq = 0.0f;
    float       fQ = 0.0f;
    if( !bLinearFit(&afSsquared[0], 
                    &afLogIOverSigmaValues[0], 
                    afLogIOverSigmaValues.size()-1, 
                    &afWeights[0], 
                    true, 
                    &fA,
                    &fB, 
                    &fSigA, 
                    &fSigB,
                    &fChiSq, 
                    &fQ) )
         return false;
        
    double      fExtrapolatedIOverSigmaAtTargetResolution = pow(10.0, (fA +
                                                            fB * (0.25/(float)fTargetReso/(float)fTargetReso) ) );
    
    Crotation           oRot(*poHeader);
    if( !oRot.bIsAvailable() ) // Rotation information not available
        return false;
    
    double      fImageRotationWidth = oRot.fGetIncrement();
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Figure out minimum required exposure time to achieve the target I/Sigma at the target resolution
    double  fImageExposureTime = oRot.fGetExposureTime(); // current exposure time
    
    double      fProposedTime = pow(fImageExposureTime, fIOverSigmaTimePower) *
                                fTargetIOverSigma/fExtrapolatedIOverSigmaAtTargetResolution;
    fProposedTime = pow(fProposedTime, 1.0/fIOverSigmaTimePower);
    fProposedTime = dGetNearestMultipleInteger(fProposedTime, 5, DTREK_VEC_ROUNDOFF_CEILING);
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Record exposure time evaluation parameters
    char    cTemp[256];
    sprintf(cTemp, "%.3f %.3f %.3f %.3f %.2f %.1f %.3f", fA,  
                                                         fSigA, 
                                                         fB, 
                                                         fSigB,
                                                         fImageRotationWidth,
                                                         fImageExposureTime,
                                                         fIOverSigmaTimePower);
    
    poHeader->nReplaceValueDictionary(strCallingModuleName, Cstring(D_K_ExposureTimeEvaluationInfo), cTemp);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Record the exposure time evaluation results
    sprintf(cTemp, "%.0f %.2f %.2f %.1f %.1f", fProposedTime, 
                                               fImageRotationWidth,
                                               fTargetReso,
                                               fTargetIOverSigma,
                                               fExtrapolatedIOverSigmaAtTargetResolution);
    poHeader->nReplaceValueDictionary(strCallingModuleName, Cstring(D_K_ProposedExposureTime), cTemp);

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void Creflnlist::vRemoveReflnsInResolutionRingAreas(Cimage_header* poHeader, double dRingWidthToExclude, DTREK_WORD wCtrl)
{
    if( !poHeader )
        return; // safety
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get ring locations from the header
    itr<double>     adResoValues;
    double          dTemp = 0.0;
    int             iRingIndex = 0;
    
    Cstring     strCallingModuleName = sDtrekGetModuleName();
    strCallingModuleName.upcase();

    strCallingModuleName += "_RANK";
    
    while( 0 == poHeader->nGetValueDictionary(strCallingModuleName, "RESO_RING_MID_RESO_*", dTemp, iRingIndex) )
    {
        adResoValues + dTemp;
        iRingIndex++;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( 0 == adResoValues.size() )
        return;  // no rings info in the header

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get wavelength from the header
    Cstring             strSourcePrefix("SOURCE_");
    Cwavelength         oWavelength(*poHeader, strSourcePrefix, false);

    double              dWavelength = (double)oWavelength.fGetWavelength();

    nDeleteRing(&adResoValues[0], 
                adResoValues.size(),
                dRingWidthToExclude,
                dWavelength,
                wCtrl);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Creflnlist::vSetColorsForDisplay()
{
    if( m_nFI_nNonunfFlag < 0 )
        m_nFI_nNonunfFlag = nExpandGetField(ms_snNonunfFlag);
    
    double      dIntensity = 0.0;
    
    int         nColor = -1;

    for(int nRef = 0; nRef < m_nNumReflns; nRef++) 
    {
        dIntensity = (*this)[nRef].fGetIntensity();
    
        if( -999.0 == dIntensity )  // -999.0 is what the intensity is set to in Crefine::nGetReflnsFromImages()
            nColor = (int)enColorRed;
        else if ( 0.0 == dIntensity )
            nColor = (int)enColorGreen;
        else
            nColor = (int)enColorBlue;
    
        (*this)[nRef].vSetField(m_nFI_nNonunfFlag, nColor);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Creflnlist::vGetResoLimits(double& fResoMin, double& fResoMax, int nResoFieldIndex)
{
    if( nResoFieldIndex < 0 )
        return; // safety
    
    double      fReso = 0.0;

    fResoMin = -1.0;
    fResoMax = -1.0;

    for(int nRef = 0; nRef < m_nNumReflns; nRef++)
    {
        if( (*this)[nRef].m_nSelect == m_pcSelect[m_nTotalFieldsPlus1-1] )
        {
            fReso = (*this)[nRef].fGetField(nResoFieldIndex);
        
            if( fResoMin == -1.0 ) 
            {
                fResoMin = fReso;
                fResoMax = fReso;
            }
            else
            {
                fResoMin = max(fResoMin,fReso);
                fResoMax = min(fResoMax,fReso);
            }
        }
    }
    
    if( fResoMin < fResoMax )
      std::swap(fResoMin, fResoMax);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Creflnlist::bTwoReflnsOnSameFrame(int index1, int index2, double fScanRotStart, double fRotWidth)
{
    double      fReflnRotStart_1 = (*this)[index1].fGetField(m_nFI_fCalcRotStart);
    double      fReflnRotEnd_1   = (*this)[index1].fGetField(m_nFI_fCalcRotEnd);
    double      fReflnRotMid_1   = (*this)[index1].fGetField(m_nFI_fCalcRotMid);

    double      fReflnRotStart_2 = (*this)[index2].fGetField(m_nFI_fCalcRotStart);
    double      fReflnRotEnd_2   = (*this)[index2].fGetField(m_nFI_fCalcRotEnd);
    double      fReflnRotMid_2   = (*this)[index2].fGetField(m_nFI_fCalcRotMid);

    double  fNumFramesFromScanStartToMidPointBetweenPeaks = 0.0;
    double	fCriticalFrameStart = 0.0;
    double	fCriticalFrameEnd = 0.0;

    if( fReflnRotMid_1 < fReflnRotMid_2 && fReflnRotEnd_1 < fReflnRotStart_2 )  
    {                                                                           
        //     *     *                                                                              
        //***** ***** ***** 
        //     1     2      
        fNumFramesFromScanStartToMidPointBetweenPeaks = ((fReflnRotEnd_1 + fReflnRotStart_2)/2.0 - fScanRotStart) / fRotWidth;
    
        fCriticalFrameStart = fScanRotStart + floor(fNumFramesFromScanStartToMidPointBetweenPeaks) * fRotWidth;
        fCriticalFrameEnd   = fScanRotStart + ceil(fNumFramesFromScanStartToMidPointBetweenPeaks) * fRotWidth;
    
        if( !(fReflnRotEnd_1 > fCriticalFrameStart && fReflnRotStart_2 < fCriticalFrameEnd) )
            return false; // reflections are NOT present on the same osc frame!
    }
    
    if( fReflnRotMid_1 > fReflnRotMid_2 && fReflnRotStart_1 > fReflnRotEnd_2 )    
    {                                                                             
        //     *     *                                                                               
        //***** ***** *****
        //     2     1     
        fNumFramesFromScanStartToMidPointBetweenPeaks = ((fReflnRotEnd_2 + fReflnRotStart_1)/2.0 - fScanRotStart) / fRotWidth;

        fCriticalFrameStart = fScanRotStart + floor(fNumFramesFromScanStartToMidPointBetweenPeaks) * fRotWidth;
        fCriticalFrameEnd   = fScanRotStart + ceil(fNumFramesFromScanStartToMidPointBetweenPeaks) * fRotWidth;

        if( !(fReflnRotEnd_2 > fCriticalFrameStart && fReflnRotStart_1 < fCriticalFrameEnd) )
            return false; // reflections are NOT present on the same osc frame!
    }
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove all reflections that have "bad" intensity or sigma value
void Creflnlist::vRemoveIntensityOrSigmaMarkedReflections()
{
    int*        pnDelete = new int[m_nNumReflns];
    
    double      fIntensity = 0.0;
    double      fSigma = 0.0;
    

    for(int ii=0; ii < m_nNumReflns; ii++)
    {
        fIntensity = (*this)[ii].fGetIntensity();
        fSigma     = (*this)[ii].fGetSigmaI();

        pnDelete[ii] = ( 0.0 >= fIntensity || 0.0 >= fSigma ) ? 1 : 0;
    }
    
    int     nSaveNumReflns = m_nNumReflns;
    
    nDelete(1, pnDelete);

    delete [] pnDelete;
    pnDelete = NULL;

    if( nSaveNumReflns > m_nNumReflns )
        printf("\n\nDeleted %d reflections not found on images(s).\n\n", nSaveNumReflns - m_nNumReflns);

}
/////////////////////////////////////////////////////////////////////////////////////////////////
// Check for overlap and remove all but one reflection from each group of overlapped reflections.
void Creflnlist::vRemoveOverlapReflns(double fRotWidth)
{
    if( fRotWidth <= 0.0 )
        return; // safety check
    
    //////////////////////////////////////////////////////////////////
    int     nFI_nHead  = nExpandGetField("nMergeHead");
    int     nFI_nNext  = nExpandGetField("nMergeNext");
        
    int     nOverlapped = 0;
    if( -1 == (nOverlapped=nOverlapCheck(nFI_nHead, nFI_nNext, fRotWidth)) )
    {
        printf("\n\nError: reflection overlap check failed.\n\n");
        return;
    }
    else
        printf("\nFound %d overlapped reflections.\n", nOverlapped);
    
    int     nSaveNumReflns = m_nNumReflns;

    int*        pnDelete = new int[m_nNumReflns];

    int     jj = 0;
    for(int ii=0; ii < m_nNumReflns; ii++)
    {
        // If this reflection is a "head" reflection, leave it, but delete all "next" reflections if any.
        // If this reflection is not a "head" reflection, skip it for now, because we will handle it when we come to its "head" reflection.
        if( ii != (*this)[ii].nGetField(nFI_nHead) )
            continue;

        pnDelete[ii] = 0;  // we are not deleting any "head" reflections
        if( -1 != (*this)[ii].nGetField(nFI_nNext) ) // a "head" reflection ii overlaps with at least one other reflection
        {
            //go through the linked list and mark all "next" reflections for deletion
            jj = ii;
            while( -1 != (jj = (*this)[jj].nGetField(nFI_nNext)) )
            {
                pnDelete[jj] = 1;
            }
        }
    }
    
    nDelete(1, pnDelete);

    delete [] pnDelete;
    pnDelete = NULL;

    if( nSaveNumReflns > m_nNumReflns )
        printf("Deleted %d overlapped reflections (kept one in each overlapped group).\n\n", nSaveNumReflns - m_nNumReflns);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Creflnlist::vToStratResoBins(CResoStrategyBins* poResoBins)
{
    if( !poResoBins )
        return;  // safety
    
    CResoStrategyBin*       pBin = NULL;
    Crefln*         poRefln = NULL;
    double          fReso = 0.0;
    int             nReject = 0;
    
    int             nFI_fReso = nGetFieldIndex("fResolution"); 
    
    if( nFI_fReso < 0 )
        return;

    for(int ii = 0; ii < nGetNumReflns(); ii++)
    {
        poRefln  = poGetRefln(ii);
        
        fReso = poRefln->fGetField(nFI_fReso);

        pBin = poResoBins->pGetBinReso(fReso);

        if( NULL == pBin ) // We need to skip that reflection
        {
            continue; 
        }

        pBin->vAddRefl();

        nReject = poRefln->nGetField(m_nFI_nNonunfFlag);

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
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
*_RESO_*_NUMSPOTS
*_IOVERSIG_AVG
*_SHARPNESS
*_FRACTION_ACCT
*/
