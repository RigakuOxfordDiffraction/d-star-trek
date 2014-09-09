#ifndef DT_DTREK_H
#define DT_DTREK_H
//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//84
// Dtrek.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the general header file for d*trek modules.
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

//+Include files

#if !defined(VC9) && !defined(__APPLE__) && !defined(NEW_IOS)
    #include <iostream.h>
    #include <iomanip.h>
	#include <fstream.h>
#else
	#include <iostream>
    #include <iomanip>
	#include <fstream>
	#include <string>
	
        // don't do this! see: http://www.parashift.com/c++-faq-lite/coding-standards.html#faq-27.5
        // using namespace std;

	//#define getline dt_getline
	//#define strncpy dt_strncpy
	//#define strcpy dt_strcpy
	//#define sscanf dt_sscanf
	//#define sprintf dt_sprintf
	//#define strtok dt_strtok
#endif

#include <map>
#include <vector>
 
#pragma warning(disable:4786)
#pragma warning(disable:4251)

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "minmax.h"

#ifdef NOANSIMATH
#include "ansimath.h"
#endif

//+Definitions and constants

// added ndef for alpha compilation since don't need macro on alpha - tlh
#if (!defined(ALPHA_VMS)) && (!defined(ICC)) && (!defined(SUNOS))
#define abs(f)			( (f) >=0 ? (f) : -(f))
#endif

#define ABS(f)			( (f) >=0 ? (f) : -(f))

#ifndef TRUE
#define TRUE          1
#endif

#ifndef FALSE
#define FALSE         0
#endif

#ifndef true
#define true          1
#endif

#ifndef false
#define false         0
#endif

#ifndef NULL
#define NULL          0
#endif

// Force support for other detector types to be included by default

#define RAXIS        1
#define SIEMENS      1
#define MARCCD       1
#define MARIP        1
#define ADSCCCD      1
#define BRANDEISCCD  1
#define DIP2030      1
#define MEDOPTICS    1

//RB Disable CBF for Windows until we sort it out
//#ifndef WIN32
    #define DTREK_CBF          1
//#endif

// Typedefs

// added typedef for bool on alpha - tlh
// also for Windows NT
#if defined(ALPHA_VMS) || defined(__NUTC__) || defined(DT_NOBOOL)
typedef unsigned char bool;
#elif defined (NEED_BOOL)
#ifndef WIN32
#include "bool.h"
#endif
#endif

#ifdef OSF1
typedef unsigned int UINT4;
typedef int           INT4;
typedef unsigned int ULONG;
typedef int           LONG;
#define atoll atol
#else
typedef unsigned int UINT4;
typedef int           INT4;
typedef unsigned long int ULONG;
typedef long int           LONG;
#endif


#ifdef THROW_INSTEAD_OF_EXIT
#define exit throw
#endif

#endif   // DT_DTREK_H

