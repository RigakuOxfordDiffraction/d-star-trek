#ifndef DT_CSPACEGROUP_H
#define DT_CSPACEGROUP_H
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
// Cspacegroup.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Cspacegroup
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

#include "Cspacegroup_table_decl.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eSpacegroup_states {
  eSpacegroup_unknown_state,
  eSpacegroup_known_state
};

enum eSpacegroup_classes {
  eSpacegroup_unknown_class,
  eSpacegroup_triclinic_class,
  eSpacegroup_monoclinic_class,
  eSpacegroup_orthorhombic_class,
  eSpacegroup_tetragonal_class,
  eSpacegroup_trigonal_class,
  eSpacegroup_hexagonal_class,
  eSpacegroup_rhombohedral_class,
  eSpacegroup_cubic_class
};

//+Code begin

class Crefln;      // Forward declaration of class Crefln
class Claue;

class DTREK_EXPORT Cspacegroup {
friend class Claue;

protected:

eSpacegroup_states m_eThe_State;       // The state of this objectn
eSpacegroup_classes m_eThe_Class;      // Class of this spacegroup
int                m_nNumber;          // Intntl tables number of the spacegroup
Cstring            m_sName;            // Intntl tables name of the spacegroup
Cstring            m_sEquiv;           // Equivalent positions
Cstring            m_sDescription;     // Description of spacegroup
int                m_nNumEquiv;        // Number of equivalent positions
int                m_nNumLaue;         // Number of laue symops, no inversion
int                m_nCentFlag;        // Centering flag 0=P, 1=A, 2=B, 3=C,
                                       //    4=F, 5=I, 6=R
float             *m_pfRotT;           // Rotation components (3x3 matrices) transposed!
float             *m_pfTrans;          // Translation components

//

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cspacegroup ();                 // Construct an empty spacegroup object

Cspacegroup (const int nNum);   // Construct spacegroup object with given number

Cspacegroup (const Cstring& sNameIn);   // Construct spacegroup object with given name

~Cspacegroup ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

int  nList(const int nFlag = 0);
bool bIsExtinct(const int nH, const int nK, const int nL);
eSpacegroup_classes eGetClass(void);
Cstring sGetClass(void);
char cGetClass(void);
int  nReadEquivPos();
int nReadEquivPos(Cstring& sDef);
Cstring sNameFromNumber(void);
int  nNumberFromName(void);
int  nNumberFromNameCaseInsensitive(void);
int  nReduceHKL(Crefln *poRefln, int *pHKLOut, int *pnPlusMinus = (int*)NULL);
void vSet(const int nNumberIn);
void vSet(const Cstring &sNameIn);
void vSetCentFlagFromName(void);
char cGetCentFlag(void);
Cstring sGetEquivPos(const int nPosn, const int nNoTransFlag=0);

inline
int  nGet(void) { return (m_nNumber);}

inline
int  nGetNumEquiv(void) { return (m_nNumEquiv); }

inline
int  nGetNumLaue(void) { return (m_nNumLaue); }

inline
int  nGetCentFlag(void) { return (m_nCentFlag); }

inline
Cstring sGetName(void)  { return (m_sName); }
//

private:

int nInitValues(void);

};  // end of class Cspacegroup

#endif   // DT_CSPACEGROUP_H

