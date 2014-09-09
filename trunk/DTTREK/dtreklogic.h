#ifndef DT_DTREKLOGIC_H
#define DT_DTREKLOGIC_H
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
// dtreklogic.h        Initial author: J.W. Pflugrath           14-Jan-1996
//    This file is the header file for implementation of some static functions
//    which perform a little boolean logic.
//
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

#include "Dtrek.h"
#include "Cstring.h"

// Static function prototypes

DTREK_EXPORT bool bEqual(const float fA, const float fB);
DTREK_EXPORT bool bNotEqual(const float fA, const float fB);
DTREK_EXPORT bool bGreaterEqual(const float fA, const float fB);
DTREK_EXPORT bool bGreater(const float fA, const float fB);
DTREK_EXPORT bool bLess(const float fA, const float fB);
DTREK_EXPORT bool bLessEqual(const float fA, const float fB);

DTREK_EXPORT bool bEqual(const int nA, const int nB);
DTREK_EXPORT bool bNotEqual(const int nA, const int nB);
DTREK_EXPORT bool bGreaterEqual(const int nA, const int nB);
DTREK_EXPORT bool bGreater(const int nA, const int nB);
DTREK_EXPORT bool bLess(const int nA, const int nB);
DTREK_EXPORT bool bLessEqual(const int nA, const int nB);
DTREK_EXPORT bool bAnd(const int nA, const int nB);
DTREK_EXPORT bool bOr(const int nA, const int nB);

DTREK_EXPORT bool bEqual(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bEqualRegular(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bNotEqual(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bNotEqualRegular(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bGreaterEqual(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bGreater(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bLess(const Cstring sA, const Cstring sB);
DTREK_EXPORT bool bLessEqual(const Cstring sA, const Cstring sB);

#endif   // DT_DTREKLOGIC_H
