#ifndef DT_CREFLN_H
#define DT_CREFLN_H
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
// Crefln.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Crefln
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

#include "Cstring.h"
#include "Creflnlist.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

static char cEOL = '\n';

#define nMaxErrorStatusTypes 32			// Status types for each possible type of error that could happen durring integration. (see ms_asErrorStrings m_anNumReflns)

enum eRefln_states {
  eRefln_unknown_state,
  eRefln_picked_state,
  eRefln_indexed_state,
  eRefln_predicted_state,
  eRefln_refined_state,
  eRefln_building_3Ddata_state,
  eRefln_ready_to_integrate,
  eRefln_integrated_state,
  eRefln_correct_state
};

//+Code begin

class Creflnlist;    // Forward reference of class Creflnlist
class Cimage_header; // Forward reference of class Cimage_header

const int nCREFLN_DATA_SIZE = 256;
const int nCREFLN_STRING_SIZE = 8;
const int nCREFLN_MAX_FIELDS = nCREFLN_DATA_SIZE/sizeof(int);

union DTREK_EXPORT CreflnData {
    int m_pnIntField[nCREFLN_DATA_SIZE/sizeof(int)];
    float m_pfFloatField[nCREFLN_DATA_SIZE/sizeof(float)];
    char m_pcStringField[nCREFLN_DATA_SIZE/nCREFLN_STRING_SIZE][nCREFLN_STRING_SIZE];
};


class DTREK_EXPORT Crefln {

friend class Creflnlist;

protected:
    char m_a2cOffsets[2];
    int  m_nTotalDataSize;

    CreflnData* m_poData;


  eRefln_states m_eThe_State;      // The enumerated State of this refln
  int nComputeOffsets(int nInts,int nFloats,int nStrings);

public:
  int           m_nSelect;
  void         *m_pvUserData;      // Pointer to user data for this reflection
  Creflnlist   *m_poReflnlist;     // Pointer to Creflnlist object to which
                                   //   this reflection belongs.

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments
////////////////////////////////////////////////////////////////////////

Crefln (Creflnlist* poReflnlistIn,
	const int nHval=0, const int nKval=0, const int nLval=0);
Crefln (const Crefln& oOther);          // Copy constructor 1
Crefln (Creflnlist& oReflnlist, const Crefln& oOther); // Copy constr. 2
Crefln& operator=(const Crefln& oOther);

~Crefln ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
////////////////////////////////////////////////////////////////////////

int  nInitValues(void);
int  nList(const int nType = 0);
int  nUpdateHeader(Cimage_header *poHeader, const Cstring& sPre="");
void vWrite(FILE* pFOut, const char *pcSelectIn = NULL);
int  nRead(FILE* pFIn);
int  nPackHKL(const int *pnHKL=NULL) const;
int  nUnPackHKL(const int nPacked, int *pnHKL=NULL);
void vReindex(const float *pfMatIn);
int nConvertToAxialEllipse(double* pfMaxSize = NULL,double* pfMinSize = NULL);
bool bEllipsoidAvailable(bool bCheckParamFieldsAlso = false);
int nPutGetEllipse(bool bPutToRefln,double* a2x2fEllipseA,double* a2fEllipseb,double* pfEllipsec,double* a2fOffset = NULL);
int nPutGetEllipse(bool bPutToRefln,float* a2x2fEllipseA,float* a2fEllipseb,float* pfEllipsec,float* a2fOffset = NULL);

inline int* pnGetHKL(void) const { return (int*) ((&m_poData->m_pnIntField[0])); }
inline int   nGetH(void) const   { return (*m_poData->m_pnIntField);       }
inline int   nGetK(void) const   { return (*(m_poData->m_pnIntField+1));   }
inline int   nGetL(void) const   { return (*(m_poData->m_pnIntField+2));   }

inline void vSetH(const int nHval)  {  *m_poData->m_pnIntField    = nHval; }
inline void vSetK(const int nKval)  { *(m_poData->m_pnIntField+1) = nKval; }
inline void vSetL(const int nLval)  { *(m_poData->m_pnIntField+2) = nLval; }

inline float fGetIntensity(void) const { return (m_poData->m_pfFloatField[m_a2cOffsets[0]]);   }
inline float fGetSigmaI(void) const    { return (m_poData->m_pfFloatField[m_a2cOffsets[0]+1]); }

inline void vSetIntensity(const float fIntval) {m_poData->m_pfFloatField[m_a2cOffsets[0]] = fIntval; }
inline void vSetSigmaI(const float fSigval)    {m_poData->m_pfFloatField[m_a2cOffsets[0]+1] = fSigval; }

inline const int nGetField(const int nFieldIndex) const
  {
    return (m_poData->m_pnIntField[nFieldIndex]);
  }

inline const float fGetField(const int nFieldIndex) const
  {
    return (m_poData->m_pfFloatField[m_a2cOffsets[0] + nFieldIndex]);
  }

inline const Cstring sGetField(const int nFieldIndex) const
  {
    return (Cstring) (&m_poData->m_pcStringField[m_a2cOffsets[1] + nFieldIndex][0]);
  }

inline const char* pcGetField(const int nFieldIndex) const
{
    return &m_poData->m_pcStringField[m_a2cOffsets[1] + nFieldIndex][0];
};

inline void vSetField(const int nFieldIndex, const int nVal)
  {
    m_poData->m_pnIntField[nFieldIndex] = nVal;
  }

inline void vSetField(const int nFieldIndex, const float fVal)
  {
    m_poData->m_pfFloatField[m_a2cOffsets[0] + nFieldIndex] = fVal;
  }
inline void vCheckStrings();


void vSetField(const int nFieldIndex, const Cstring& sVal);
void vSetField(const int nFieldIndex, const char* pcVal);

  static int ms_nErrorSpecial;
  static int ms_nErrorOnEdge0;
  static int ms_nErrorOffEdge0;
  static int ms_nErrorOnEdge1;
  static int ms_nErrorOffEdge1;
  static int ms_nErrorOnEdge2;
  static int ms_nErrorOffEdge2;
  static int ms_nErrorHKLBound;
  static int ms_nErrorRing;
  static int ms_nErrorTooDark;
  static int ms_nErrorBackground;
  static int ms_nErrorBackgroundSD;
  static int ms_nErrorNonunfA;
  static int ms_nErrorNonunfB;
  static int ms_nErrorPartialStart;
  static int ms_nErrorPartialEnd;
  static int ms_nErrorRotTooWide;
  static int ms_nErrorIntProblem;
  static int ms_nErrorTooFewPix;
  static int ms_nErrorUnaccountedPix;
  static int ms_nErrorOverlap;

DTREK_WIN_DLL_DATA_EXPORT  static Cstring ms_asErrorsNames[nMaxErrorStatusTypes];    // Command line names of errors used mainly with '-rogue' option.
  static Cstring ms_asErrors[nMaxErrorStatusTypes];         // nMaxErrorStatusTypes possible error strings
  static Cstring ms_asErrorPhrases[nMaxErrorStatusTypes];   // nMaxErrorStatusTypes possible error phrases
  static Cstring ms_asErrorStrings[nMaxErrorStatusTypes];   // Error messages

};  // end of class Crefln

typedef struct
{
	void vInitialize()
	{
		m_nSaturatedCount=0;			
		m_dMaxPixelValue=0.0;			
		m_nPeakEllipsoidPixelCount=0;	
	}
	int			m_nSaturatedCount;						// number of saturated pixels for peak
	double		m_dMaxPixelValue;						// maximum pixel value for peak
	int			m_nPeakEllipsoidPixelCount;				// number of pixels in the peak elliptical mask
}REFLN_EXTRA_INFO;

#endif   // DT_CREFLN_H

