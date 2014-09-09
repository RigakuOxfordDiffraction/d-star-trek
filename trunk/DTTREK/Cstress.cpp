//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cstress.cpp       Initial author: RB     09-Nov-2006
// This file contains the member functions of class Cstress

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
 
#include "Cstress.h"

Cstress::Cstress()
{
    vInit();
}

Cstress::Cstress(Cimage_header& oHeader)
{
    vInit();

    // Read number of Z positions
    oHeader.nGetValue(Cstring(D_K_StressNumOfZPositions), &m_nNumOfZPositions);
    oHeader.nGetValue(Cstring(D_K_StressSlitWidth), &m_nSlitWidth);
    oHeader.nGetValue(Cstring(D_K_StressInclinationType), &m_sInclinationType);
    
    oHeader.nGetValue(Cstring(D_K_StressZPositions), m_anZPositions, " ");
    oHeader.nGetValue(Cstring(D_K_StressPsiAngles), m_afPsiAngles, " ");

    oHeader.nGetValue(Cstring(D_K_StressPeakAngle), &m_fPeakAngle);
    oHeader.nGetValue(Cstring(D_K_StressYoungModulus), &m_fYoungModulus);
    oHeader.nGetValue(Cstring(D_K_StressPoissonRatio), &m_fPoissonRatio);
    oHeader.nGetValue(Cstring(D_K_StressStressConstant), &m_fStressConstant);
}

Cstress::Cstress(const Cstress& oOther)
{
    m_nNumOfZPositions = oOther.m_nNumOfZPositions;
    m_nSlitWidth       = oOther.m_nSlitWidth;       
    m_sInclinationType = oOther.m_sInclinationType;
  
    int     ii = 0;
    m_anZPositions.clear();
    for(ii=0; ii < oOther.m_anZPositions.size(); ii++)
        m_anZPositions.push_back(oOther.m_anZPositions[ii]);

    m_afPsiAngles.clear();
    for(ii=0; ii < oOther.m_afPsiAngles.size(); ii++)
        m_afPsiAngles.push_back(oOther.m_afPsiAngles[ii]);

    m_fPeakAngle        = oOther.m_fPeakAngle;	        
    m_fYoungModulus     = oOther.m_fYoungModulus;         
    m_fPoissonRatio     = oOther.m_fPoissonRatio;         
    m_fStressConstant   = oOther.m_fStressConstant;       
}

Cstress& Cstress::operator=(const Cstress& oOther)
{
    m_nNumOfZPositions = oOther.m_nNumOfZPositions;
    m_nSlitWidth       = oOther.m_nSlitWidth;       
    m_sInclinationType = oOther.m_sInclinationType;
  
    int     ii = 0;
    m_anZPositions.clear();
    for(ii=0; ii < oOther.m_anZPositions.size(); ii++)
        m_anZPositions.push_back(oOther.m_anZPositions[ii]);

    m_afPsiAngles.clear();
    for(ii=0; ii < oOther.m_afPsiAngles.size(); ii++)
        m_afPsiAngles.push_back(oOther.m_afPsiAngles[ii]);

    m_fPeakAngle        = oOther.m_fPeakAngle;	        
    m_fYoungModulus     = oOther.m_fYoungModulus;         
    m_fPoissonRatio     = oOther.m_fPoissonRatio;         
    m_fStressConstant   = oOther.m_fStressConstant;       

    return *this;
}

Cstress::~Cstress()
{
}

void Cstress::vInit()
{
    m_nNumOfZPositions = 0;
    m_nSlitWidth       = 0;       
    m_sInclinationType = "side";

    m_anZPositions.clear();
    m_afPsiAngles.clear();
  
    m_fPeakAngle        = 0.0;	        
    m_fYoungModulus     = 0.0;         
    m_fPoissonRatio     = 0.0;         
    m_fStressConstant   = 0.0;       
}

int Cstress::nUpdateHeader(Cimage_header& oHeader)
{
    oHeader.nReplaceValue(Cstring(D_K_StressNumOfZPositions), m_nNumOfZPositions);
    oHeader.nReplaceValue(Cstring(D_K_StressSlitWidth), m_nSlitWidth);
    oHeader.nReplaceValue(Cstring(D_K_StressInclinationType), m_sInclinationType);
    
    if( m_anZPositions.size() > 0 )
        oHeader.nReplaceValue(Cstring(D_K_StressZPositions), m_anZPositions.size(), &m_anZPositions[0]);
    else
        oHeader.nReplaceValue(Cstring(D_K_StressZPositions), Cstring(""));

    if( m_afPsiAngles.size() > 0 )
        oHeader.nReplaceValue(Cstring(D_K_StressPsiAngles), m_afPsiAngles.size(), &m_afPsiAngles[0]);
    else
        oHeader.nReplaceValue(Cstring(D_K_StressPsiAngles), Cstring(""));

    oHeader.nReplaceValue(Cstring(D_K_StressPeakAngle),      m_fPeakAngle);
    oHeader.nReplaceValue(Cstring(D_K_StressYoungModulus),   m_fYoungModulus);
    oHeader.nReplaceValue(Cstring(D_K_StressPoissonRatio),   m_fPoissonRatio);
    oHeader.nReplaceValue(Cstring(D_K_StressStressConstant), m_fStressConstant);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////
