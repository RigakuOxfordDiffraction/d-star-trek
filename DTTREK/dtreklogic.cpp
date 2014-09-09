//
// Copyright (c) 1996 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtreklogic.cc     Initial author: J.W. Pflugrath           13-Jan-1996
//    This file contains some logical static functions
//    that do not belong to any class.
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

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "dtreklogic.h"           // Function prototypes

//+Definitions, constants, and initialization of static member variables

//+Code begin

//+Public functions

bool bEqual(const float fA, const float fB)
{ return (fA == fB); }

bool bNotEqual(const float fA, const float fB)
{ return (fA != fB); }

bool bGreaterEqual(const float fA, const float fB)
{ return (fA >= fB); }

bool bGreater(const float fA, const float fB)
{ return (fA > fB); }

bool bLess(const float fA, const float fB)
{ return (fA < fB); }

bool bLessEqual(const float fA, const float fB)
{ return (fA <= fB); }


bool bEqual(const int nA, const int nB)
{ return (nA == nB); }

bool bNotEqual(const int nA, const int nB)
{ return (nA != nB); }

bool bGreaterEqual(const int nA, const int nB)
{ return (nA >= nB); }

bool bGreater(const int nA, const int nB)
{ return (nA > nB); }

bool bLess(const int nA, const int nB)
{ return (nA < nB); }

bool bLessEqual(const int nA, const int nB)
{ return (nA <= nB); }

bool bAnd(const int nA,const int nB)
{ return ((((unsigned int) nA) & ((unsigned int) nB)) != 0); };

bool bOr(const int nA,const int nB)
{ return ((((unsigned int) nA) | ((unsigned int) nB)) != 0); };


bool bEqual(const Cstring sA, const Cstring sB)
{ return (sA == sB); }

bool bNotEqual(const Cstring sA, const Cstring sB)
{ return (sA != sB); }

bool bEqualRegular(const Cstring sA, const Cstring sB) 
{ return (sA.bRegularMatch(sB)); };

bool bNotEqualRegular(const Cstring sA,const Cstring sB)
{ return (!sA.bRegularMatch(sB)); };

bool bGreaterEqual(const Cstring sA, const Cstring sB)
{ return (sA >= sB); }

bool bGreater(const Cstring sA, const Cstring sB)
{ return (sA > sB); }

bool bLess(const Cstring sA, const Cstring sB)
{ return (sA < sB); }

bool bLessEqual(const Cstring sA, const Cstring sB)
{ return (sA <= sB); }
