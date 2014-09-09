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
// ansimath.cc     Initial author: J.W. Pflugrath           01-Sep-1995
//    This file contains some math float functions not found on
//    some platforms.
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
// Some ANSI math items that may be missing on some systems

#include <math.h>
#include "ansimath.h"

float sqrtf(const float fTemp)
{
  return ( (float) sqrt((double) fTemp));
}
/*************************
The following should NOT be used by d*TREK modules

float cosf(const float fTemp)
{
  return ( (float) cos((double) fTemp));
}

float sinf(const float fTemp)
{
  return ( (float) sin((double) fTemp));
}

float tanf(const float fTemp)
{
  return ( (float) tan((double) fTemp));
}

float acosf(const float fTemp)
{
  return ( (float) acos((double) fTemp));
}

float asinf(const float fTemp)
{
  return ( (float) asin((double) fTemp));
}

float atan2f(const float fTemp1, const float fTemp2)
{
  return ( (float) atan2((double) fTemp1, (double) fTemp2));
}

float atanf(const float fTemp)
{
  return ( (float) atan((double) fTemp));
}

float powf(const float fTemp1, const float fTemp2)
{
  return ( (float) pow((double) fTemp1, (double) fTemp2));
}

float logf(const float fTemp)
{
  return ( (float) log((double) fTemp));
}

float log10f(const float fTemp)
{
   return( (float) log10((double) fTemp));
}

float fabsf(const float fTemp)
{
   return( (float) fabs((double) fTemp));
}
****************************/
