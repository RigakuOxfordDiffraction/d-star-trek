//
// Copyright (c) 2006 Rigaku
//
//
// Cstress.h        Initial author: S.Yasukawa           03-Nov-2006
// This file is the header file for class Cstress

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
 
#ifndef DT_CSTRESS_H
#define DT_CSTRESS_H

//+Include files

#include "Cimage_header.h"

//+Definitions and constants


//+Code begin

class DTREK_EXPORT Cstress
{
public:
protected:
  int    m_nNumOfZPositions;        // Number of Z positions from the header. NOTE: the m_anZPositions 
                                    // array may have a different number of Z positions.
  int    m_nSlitWidth;              // Slit width to processing on the image data
  
  Cstring m_sInclinationType;       // Inclination type. It should be "Side" or "Iso".
  
  std::vector<int>          m_anZPositions;
  std::vector<double>       m_afPsiAngles;

  double   m_fPeakAngle;	        // Peak 2theta angle to process. 
  double   m_fYoungModulus;         // Young mudulus
  double   m_fPoissonRatio;         // Posisson ratio
  double   m_fStressConstant;       // Stress constant
    
  //bool     m_bValid;
////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments
public:

Cstress();                          // Construct an empty stress object

// Construct stress object from header :
Cstress(Cimage_header& oHeader);

// Copy constructor
Cstress(const Cstress& oStress);

~Cstress();

Cstress& operator=(const Cstress& oOther);


//////////////////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
//
//////////////////////////////////////////////////////////////////////////////////////
inline int  nGetNumOfZPositions(void) { return (m_nNumOfZPositions); } // from header
inline int nGetNumOfZPositionsInArray(){return m_anZPositions.size();}
inline int nGetZPosition(int index)
{
    if(index > m_anZPositions.size() - 1 || index < 0) 
        return 0;

    return m_anZPositions[index];
}

inline int  nGetSlitWidth(void) { return (m_nSlitWidth); }

inline Cstring sGetInclinationType(void) { return (m_sInclinationType); }
inline bool bIsSideInclination(void)
{
    return (0 == m_sInclinationType.compare_no_case("side"));    
}

int nGetNumOfPsiAnglesInArray(void){return (int)m_afPsiAngles.size();}

double fGetPsiAngle(int index)
{
    if( index > (int)m_afPsiAngles.size() - 1 || index < 0 ) 
        return 0;

    return m_afPsiAngles[index];
}

inline double fGetPeakAngle(void) { return (m_fPeakAngle); }
inline double fGetYoungModulus(void) { return (m_fYoungModulus); }
inline double fGetPoissonRatio(void) { return (m_fPoissonRatio); }
inline double fGetStressConstant(void) { return (m_fStressConstant); }
///////////////////////////////////////////////////////////////////////////////////////
inline void vSetZPositions(const int nNumOfZPositions, const int* pnZpositions)
{
    m_anZPositions.clear();
    
    m_nNumOfZPositions = nNumOfZPositions;
    
    for(int ii=0; ii < nNumOfZPositions; ii++)
        m_anZPositions.push_back(pnZpositions[ii]);
}

inline void vSetNumOfZpositions(const int nPos){ m_nNumOfZPositions = nPos; }
inline void vSetSlitWidth(const int nWidth) { m_nSlitWidth = nWidth; }
inline void vSetInclinationType(const Cstring& sInclinationType){ m_sInclinationType = sInclinationType; }

inline void vSetPsiAngles(int nNumPsiAngles, double* pfPsiAngles)
{
    m_afPsiAngles.clear();
    
    for(int ii=0; ii < nNumPsiAngles; ii++)
        m_afPsiAngles.push_back(pfPsiAngles[ii]);
}

inline void vSetPeakAngle(const double fPeakAngle){ m_fPeakAngle = fPeakAngle; }
inline void vSetYoungModulus(const double fYoungModulus){ m_fYoungModulus = fYoungModulus; }
inline void vSetPoissonRatio(const double fPoissonRatio){ m_fPoissonRatio = fPoissonRatio; }
inline void vSetStressConstant(const double fStressConstant){ m_fStressConstant = fStressConstant; }

int nUpdateHeader(Cimage_header& oHeader);

private:
    void vInit();
};  // end of class Cstress

#endif   // DT_CSTRESS_H

