//
// Copyright 1999 Molecular Structure Corporation
//                9009 New Trails Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// Cstat.h    Initial author: T.L.Hendrixson           Dec 1999
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







#ifndef Cstat_H
#define Cstat_H

#include "dtrekdefs.h"


// 
//Copyright 1999 Molecular Structure Corporation  
//               9009 New Trails Drive  
//               The Woodlands, TX, USA  77381  
// 
//The contents are unpublished proprietary source  
//code of Molecular Structure Corporation  
// 
//All rights reserved   
// 
//Cstat.h      Initial author: T.L.Hendrixson           Dec 1999
/*
 * RCS stuff:
 *   $ $
 *   $ $
 *   $ $
 *   $ $
*/

typedef double CstatValue;

typedef int CstatMarker;


// Description of class goes here


class DTREK_EXPORT Cstat
{
    static CstatValue ms_fErrRef;   // Default reference if CstatValue& is required.
    int    m_nDefAlloc;             // Default allocation size.

public:

    enum eMarkTest { S_ABOVE,S_BELOW,S_ABOVEBELOW,S_ABOVEBELOW_ITER,V_ABOVE,V_BELOW,ALL,P_ABOVE,P_BELOW };


   // Create an instance of the class.
   Cstat(int nInitialSize=0,int nDefAlloc=16);  

   Cstat(const Cstat &oCopy,bool bCopyMarkedOnly= false,int nDefAlloc=16);

   virtual ~Cstat();

   // Copy an object object

   Cstat& roCopy(const Cstat& oCopy,bool bCopyMarkedOnly = false);

   // Adds a new entry to the list.  Is marked by default.
   void vAdd(CstatValue d,CstatValue dWeight=1.0,bool bMarked=true);
   // Clear the list of entries
   void vClear(void);
   int  nSize() const { return m_nNumElements; };
   int  nSetSize(int nNewSize) { vGrowList(nNewSize); m_nNumElements=nNewSize; return 0; };
   
   // Get a value, 
   CstatValue& fGet(int nIndex) const { return ((nIndex>=m_nNumElements)?(ms_fErrRef):m_fValues[nIndex]); };
   // Is the index marked? 
   bool bIsMarked(int nIndex) const { return ((nIndex>=m_nNumElements)?(false):(m_cMarkers[nIndex]!=0));  };
   // Mark/Unmark facility.
   bool bMark(bool bMark,int nIndex=-1);

   // Statistical functions.

   // Returns the average of the elements in the list.
   CstatValue fAverage(bool bRecompute=true);

   // Returns the standard deviation of the elements in the list.
   CstatValue fStandardDeviation(bool bRecompute=true);
   CstatValue fSigma(bool bRecompute=true) { return fStandardDeviation(bRecompute); };

   // Minimum value.
   CstatValue fMin();
   // Maximum value.
   CstatValue fMax();
   // Minimum value index.
   int  nMin();
   // Maximum value index.
   int  nMax();

   // Covariance
   CstatValue fCovariance(Cstat& oOther);

   // Non-accumulating vAdd(), average and statistical deviation.
   void vAddRaw(CstatValue d,CstatValue fWeight = 1.0);
   void vClearRaw();
   CstatValue fAverageRaw();
   CstatValue fStandardDeviationRaw();
   CstatValue fSigmaRaw() { return fStandardDeviationRaw(); };




   // Deviates and standard functions.
   CstatValue fNormalDeviate(CstatValue fAverage,CstatValue fSigma);
   CstatValue fPoissonDeviate(CstatValue fAverage);
   CstatValue fUniformDeviate(int nSeed = 0);
   CstatValue fGamma(CstatValue fX);

   int nMarkStat(bool bMark,eMarkTest eTest,CstatValue fVal=0.0);

   // Operators.

   CstatValue& operator [] (int nIndex) { return fGet(nIndex); };
   Cstat& operator + (const CstatValue& fAdd) { return (vAdd(fAdd),(*this)); };
   Cstat& operator = (const Cstat& oCopy) { return roCopy(oCopy); };


protected:
private:


   // This routine increase the size of the list
   void vGrowList(int nSizeDesired);


   int m_nNumElements;      // Number of elements in list
   int m_nNumAllocated;     // Number of list element allocated
   

   CstatValue* m_fValues;    // the list
   CstatValue* m_fWeights;   // Weightings.
   CstatMarker* m_cMarkers;  // marker values;

   // These store the last values of the standard deviation and the average.
   CstatValue m_fLastAverage;
   CstatValue m_fLastDeviation;
   CstatValue m_fSumX;  // Sum X
   CstatValue m_fSumX2; // Sum X^2
   CstatValue m_fSumXW; // Weighted Sum X.
   CstatValue m_fSumW;  // Sum Weights

};

#endif /* TIMINGINFO_H */
