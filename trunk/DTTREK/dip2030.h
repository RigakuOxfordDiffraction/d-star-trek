#ifndef DT_DIP2030_H
#define DT_DIP2030_H
//
// Copyright (c) 1997 Molecular Structure Corporation
//
// dip2030.h        Initial author: J.W. Pflugrath           12-Feb-1998
//    This file is the header file for a MARCCD style image file
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

#include "Cstring.h"

//+Code begin

//  The following header definition comes from M. White, UTMB

/********************TAILER STRUCTURE**************************************/


typedef struct _tagDIP2030_CO_TAIL{
        char    id[4];              /* data mark (general = DIP0)         */
        INT4    m_type1;            /* machine reading mode:              */
                                    /* == 0:line                          */
                                    /* else:spiral(IP diameter[0.1mm]     */
        INT4    m_type2;            /* 0: flat cassette                   */
                                    /* 1: Cylindrical cassette (Para)     */
                                    /* 2: Cylindrical cassette (Perp)     */
        INT4    dtype;              /* 0: short integer(4+12bits)         */
                                    /* 1: unsigned short(16bits)          */
                                    /* 2: short integer(1+15bits)         */
                                    /* n ... ... ... ... ... ... ...      */

        INT4    pixelsize;          /* pixel size (main) (um)             */
        INT4    pixelsize2;         /* pixel size (sub)  (um)             */
        float   radius;             /* radius of cylinder (mm)            */
        INT4    xsize;              /* Number of pixels along x coord     */
        INT4    ysize;              /* Number of pixels along y coord     */
        INT4    ipno;               /* imaging plate No.(1/2)             */
        char    comment[80];        /* comment message                    */
        float   x_lamda;            /* X-ray wave length  (A)             */
        float   cdist;              /* camera distance (mm)               */
        char    monochro[32];       /*  monochro parameter                */
        float   pttheta;            /* 2Theta angle (deg)                 */
                                    /* in Weissenberg configuration       */
                                    /* this is used to define mu angle    */
        INT4    ipx;                /* P position (x) of pixel adress     */
        INT4    ipy;                /* P position (y) of pixel adress     */
        float   exposure;           /* Exposure time (sec)                */
        float   kv;                 /* X.G. voltage (kv)                  */
        float   ma;                 /* X.G. current (mA)                  */
        float   collimator;         /* collimator diameter (mm)           */
        float   coupling;           /* ==  0: No weissenberg motion       */
                                    /* not 0: Weissenberg motion          */
                                    /*    for DIP2000 unit: deg/deg       */
                                    /*    for DIP3000 orcylindrical       */
                                    /*           type unit: mm/deg        */
        float    phi1;              /* Phi start angle (deg)              */
        float    phi2;              /* Phi ended angle (deg)              */
        float    phispeed;          /* Phi speed (deg/min)                */
        INT4     repet;             /* repetition number                  */
        INT4     osc_axis;          /* oscillation axis                   */
                                    /* 0: phi                             */
                                    /* 1: omega                           */
                                    /* 2: kappa                           */
        float    g_omega;           /* omega angle (deg)                  */
        float    g_kappa;           /* kappa angle (deg)                  */
        float    g_phi;             /* phi   angle (deg)                  */
        INT4     xstart;            /* start pixel position along x coord */
        INT4     ystart;            /* start pixel position along y coord */
        char     dummy[152];        /* for fist part (384 bytes)          */

        } tagDIP2030_CO_TAIL;

typedef struct _tagDIP2030_SR_TAIL {
        INT4 colour;                    /* Current color mode             */
        INT4 max_val;                   /* Current Max_Show_Value         */
        INT4 min_val;                   /* Current Min_Show_Value         */
        INT4 start_x;                   /* AOI start point ipx value      */
        INT4 start_y;                   /* AOI start point ipy value      */
        INT4 AOI_w;                     /* AOI width (in IP pixel)        */
        INT4 AOI_h;                     /* AOI width (in IP pixel)        */
        INT4 f1x;                       /* Fiducial point 1 in x          */
        INT4 f1y;                       /* Fiducial point 1 in y          */
        INT4 f2x;                       /* Fiducial point 2 in x          */
        INT4 f2y;                       /* Fiducial point 2 in y          */
        char dummy[596];                /* part 2 for screen parameters   */
        } tagDIP2030_SR_TAIL;

typedef struct _tagDIP2030_header {
        tagDIP2030_CO_TAIL tPart1;      /* part1  284 bytes               */
        tagDIP2030_SR_TAIL tPart2;      /* part2  640 bytes               */
                                        /* Total 1024 bytes               */
} tagDIP2030_header;

// Forward declaration of classes needed by the function prototypes below

class Cimage_header;

int nReadDIP2030Header(const int* pnFile, const long lFileSize, 
		       const Cstring& rsFilename,
		       char *pcBuffer, Cimage_header *poHeader);
int nReadDIP2030Trailer(const int nFile, Cimage_header *poHeader);
#endif   // DT_DIP2030_H
