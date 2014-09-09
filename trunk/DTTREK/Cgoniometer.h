#ifndef DT_CGONIOMETER_H
#define DT_CGONIOMETER_H
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
// Cgoniometer.h        Initial author: J.W. Pflugrath      01-May-1995
//    This file is the header file for class Cgoniometer
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
#include "Cstring.h"
#include "Cimage_header.h"
#include "dtrekvec.h"
#include "Segment.h"


//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eGoniometer_states {
    eGoniometer_unknown_state,
        eGoniometer_moving_state,
        eGoniometer_not_ready_to_move_state,
        eGoniometer_ready_to_move_state,
        eGoniometer_jammed_state,
        eGoniometer_collided_state
};

enum eGoniometer_types {
    eGoniometer_unknown_type,
        eGoniometer_single_type,
        eGoniometer_eulerian_type,
        eGoniometer_kappa_type,
        eGoniometer_other_type
};

const int nMaxGonioRotationOffsets = 10;   // Maximum number of distinct goniometer offset entries allowed.

//+Code begin

class DTREK_EXPORT Cgoniometer {
    
    
public:
    
    eGoniometer_states m_eThe_State;       // The state of the goniometer
    eGoniometer_types  m_eThe_Type;        // The type of the goniometer
    int                m_nNumValues;       // Number of axes or translations of goniometer
    int                m_nNumRotValues;    // Number of axes only.
    int                m_nNumTransValues;  // Number of translations only.  m_nNumRotValues + m_nNumTransValues == m_nNumValues.
    
    float             *m_pfCurrentValue;   // Current values of axes
    
    float*             m_pfHardwareMinLimit;
    float*             m_pfHardwareMaxLimit;
    
    float             *m_pfRequestedValue;
    
    float             *m_pfDatumValue;
    float             *m_pfDatumCollisionOffset;
    
    float             *m_pfSetValue;
    Cstring           *m_psName;
    Cstring           *m_psUnits;
    float             *m_pfRequestedRate;
    float             *m_pfSetRate;
    float             *m_pfVector;
    Cstring            m_sKey;
    Cstring            m_sDescription;
    Cstring            m_sPrefix;
    Cstring            m_sRelPrefix;
    
    static Cstring     ms_sGonioNumValues;
    static Cstring     ms_sGonioKey;
    static Cstring     ms_sGonioNames;
    static Cstring     ms_sGonioUnits;
    static Cstring     ms_sGonioVectors;
    static Cstring     ms_sGonioValues;
    
    static Cstring     ms_sGonioValuesMin;
    static Cstring     ms_sGonioValuesMax;

    static Cstring     ms_sGonioCollisionOffsets;
    static Cstring     ms_sGonioDescription;
    
    
    ////////////////////////////////////////////////////////////////////////
    //+Constructors, destructors and assignments
    
    Cgoniometer ();                       // Construct an empty goniometer object
    Cgoniometer(const int nNumValues, 
        const float *pfAxes, const float *pfValues,
        const Cstring *psNames, const Cstring *psUnits);
    
    // Construct goniometer object from image header:
    
    Cgoniometer (Cimage_header& oHeader, const Cstring& sPre);
    
    ~Cgoniometer ();
    
    ////////////////////////////////////////////////////////////////////////
    //+Member Function prototypes
public:
    
    int nList(const int nFlag=0);
    inline bool bIsAvailable(void) {return (m_eThe_State != eGoniometer_unknown_state); }
    inline int nGetNumValues(void) { return (m_nNumValues); }
    
    int nGetDatum(const int nNum, float *pfDatum);
    int nSetDatum(const int nNum, const float *pfDatum);
    int nGetDatum(const int nNum, double *pdDatum);
    int nSetDatum(const int nNum, const double *pdDatum);
    
    int  nGetNames(const int nNum, Cstring *psNames);
    Cstring sGetName(const int nNum);   // Convert number to name
    int  nGetNum(const Cstring sName);  // Convert name to number
    int  nGetVectors(const int nNum, float *pfVec);
    int  nGetVectors(const int nNum, double *pdVec);
    int  nGetRotVector(const int nVecNum, float *pfVec);
    int  nGetRotVector(Cstring& sName,float *pfVec);
    void vCalcGetRotMatrix(float *pfMatrix,int nRotAxis,
        const float fRot1 = -999.0,
        const float fRot2 = -999.0,
        const float fRot3 = -999.0,
        const float fRot4 = -999.0);
    
    void vCalcGetRotMatrix(double *pfMatrix,int nRotAxis,
        const double fRot1 = -999.0,
        const double fRot2 = -999.0,
        const double fRot3 = -999.0,
        const double fRot4 = -999.0);
    
    int  nCalcGetTransVector(float *pfVector);
    int  nCalcGetTransVectorDeriv(int nTransComp,double* pfVecOut);
    int  nSetTransVector(const int nDim, const float *pfVector);
    int  nUpdateHeader(Cimage_header* poHeader, const Cstring& sPre);
    
    float fGetDistance(void);

    bool bGetHardwareLimits(const char* ccAxisName, float& fMin, float& fMax);
    
    bool bGetHardwareLimits(const char* ccAxisName, CSegment& segLimits);
    
    bool bGetCollisionOffset(const char* ccAxisName, float& fOffset);

    float fGetSwing(void);
    int  nSetDistance(const float fDist);
    int  nSetSwing(const float fSwing);
    int  nDiff(Cgoniometer& oGoniometerAdd,Cgoniometer& oGoniometerSubtract);
    int  nReconcileHeaderWithImageGoniometer(Cgoniometer* poGoniometer);
    
    
    // The XXXOffset functions work on the goniometer miss-settings.
    // The XXXRotOffset functions work on the rotation offset (as for example adjusted by nGetMosaicity())

    int nFindGonioOffsets(float* pfDatumValue,int nRotVector,int& nSet,int& nPoint);
    int nApplyAllGonioOffsets(float* pfDatumValue,int nRotVector);
    int nSaveGonioOffsets(float* pfDatumValue,int nRotVector,int nNumDatum,double* pfRotOffsets);
    int nGetGonioRotOffset(float* pfDatumValue,int nRotVector,double& fRot);
    int nSaveGonioRotOffset(float* pfDatumValue,int nRotVector,int nNumDatum,double fRot);
    int nSaveOffsets(Cimage_header* poHeader);
    int nLoadOffsets(Cimage_header* poHeader);    

    
private:
    
    double  m_afGonioOffsetTableOffset[nMaxGonioRotationOffsets*5];      // Table containing goniometer values as written/read from the header file.
    double  m_afGonioOffsetTableGonio[nMaxGonioRotationOffsets*5];       // Table containing goniometer offsets as written/read from the header file.
    int     m_anGonioOffsetTableEntries[nMaxGonioRotationOffsets+1];     // Table containing number of goniometer offset entries as written/read from the header file.
    int     m_anGonioOffsetTableRotAxis[nMaxGonioRotationOffsets];       // Rotation axis used (1==first axis).
        
    int nInitValues(void);
    int nInitValues(Cimage_header& oHeader, const Cstring& sPre);
    void vGetNumRotTransValues();
    
};  // end of class Cgoniometer

#endif   // DT_CGONIOMETER_H
