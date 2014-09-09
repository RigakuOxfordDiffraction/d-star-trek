#ifndef MSC_MINMAX_H
#define MSC_MINMAX_H

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
// minmax.h        Initial author: T.L. Hendrixson          08-Feb-1996
//    This file simply defines min() and max() functions.
// There may be comparable library functions on DEC ALPHA, but I was not
// able to find them again ( I did find them once, but I don't remember
// where ).
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

// Microsoft defines __min and __max macros which screw up our min and
// max typechecking functions, so just define then here as the macros

#ifdef WIN32

#ifndef min
#define min __min
#endif

#ifndef max
#define max __max
#endif

#else

// Minimum functions

inline char   min(char   a, char   b){return (a < b) ? a : b;}

inline short  min(short  a, short  b){return (a < b) ? a : b;}

inline int    min(int    a, int    b){return (a < b) ? a : b;}
inline int    min(short  a, int    b){return (a < b) ? (int)a : b;}
inline int    min(int    a, short  b){return (a < b) ? a : (int)b;}

inline long   min(long   a, long   b){return (a < b) ? a : b;}
inline long   min(short  a, long   b){return (a < b) ? (long)a : b;}
inline long   min(int    a, long   b){return (a < b) ? (long)a : b;}
inline long   min(long   a, short  b){return (a < b) ? a : (long)b;}
inline long   min(long   a, int    b){return (a < b) ? a : (long)b;}

inline float  min(float  a, float  b){return (a < b) ? a : b;}
inline float  min(short  a, float  b){return (a < b) ? (float)a : b;}
inline float  min(int    a, float  b){return (a < b) ? (float)a : b;}
inline float  min(long   a, float  b){return (a < b) ? (float)a : b;}
inline float  min(float  a, short  b){return (a < b) ? a : (float)b;}
inline float  min(float  a, int    b){return (a < b) ? a : (float)b;}
inline float  min(float  a, long   b){return (a < b) ? a : (float)b;}

inline double min(double a, double b){return (a < b) ? a : b;}
inline double min(short  a, double b){return (a < b) ? (double)a : b;}
inline double min(int    a, double b){return (a < b) ? (double)a : b;}
inline double min(long   a, double b){return (a < b) ? (double)a : b;}
inline double min(float  a, double b){return (a < b) ? (double)a : b;}
inline double min(double a, short  b){return (a < b) ? a : (double)b;}
inline double min(double a, int    b){return (a < b) ? a : (double)b;}
inline double min(double a, long   b){return (a < b) ? a : (double)b;}
inline double min(double a, float  b){return (a < b) ? a : (double)b;}

// Maximum functions

inline char   max(char   a, char   b){return (a > b) ? a : b;}

inline short  max(short  a, short  b){return (a > b) ? a : b;}

inline int    max(int    a, int    b){return (a > b) ? a : b;}
inline int    max(short  a, int    b){return (a > b) ? (int)a : b;}
inline int    max(int    a, short  b){return (a > b) ? a : (int)b;}

inline long   max(long   a, long   b){return (a > b) ? a : b;}
inline long   max(short  a, long   b){return (a > b) ? (long)a : b;}
inline long   max(int    a, long   b){return (a > b) ? (long)a : b;}
inline long   max(long   a, short  b){return (a > b) ? a : (long)b;}
inline long   max(long   a, int    b){return (a > b) ? a : (long)b;}

inline float  max(float  a, float  b){return (a > b) ? a : b;}
inline float  max(short  a, float  b){return (a > b) ? (float)a : b;}
inline float  max(int    a, float  b){return (a > b) ? (float)a : b;}
inline float  max(long   a, float  b){return (a > b) ? (float)a : b;}
inline float  max(float  a, short  b){return (a > b) ? a : (float)b;}
inline float  max(float  a, int    b){return (a > b) ? a : (float)b;}
inline float  max(float  a, long   b){return (a > b) ? a : (float)b;}

inline double max(double a, double b){return (a > b) ? a : b;}
inline double max(short  a, double b){return (a > b) ? (double)a : b;}
inline double max(int    a, double b){return (a > b) ? (double)a : b;}
inline double max(long   a, double b){return (a > b) ? (double)a : b;}
inline double max(float  a, double b){return (a > b) ? (double)a : b;}
inline double max(double a, short  b){return (a > b) ? a : (double)b;}
inline double max(double a, int    b){return (a > b) ? a : (double)b;}
inline double max(double a, long   b){return (a > b) ? a : (double)b;}
inline double max(double a, float  b){return (a > b) ? a : (double)b;}

#endif /* ifdef WIN32 */

#endif   /*  MSC_MINMAX_H  */
