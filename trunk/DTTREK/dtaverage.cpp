//
// Copyright (c) 1996-2006 Rigaku Americas Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtaverage.cc     Initial author: J.W. Pflugrath           07-May-1996
// Average a series of images on a pixel by pixel basis.  Exclude
// outliers from the averaging process.  Outliers are specified by the
// naming of an average image, a standard deviation image and a sigma
// value.
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
//
//    Example:  % dtaverage [ -scan headerfile] [options ...]
//                  
//    Command line options:               Default:
//                -scan     sScanFile     none
//                -template sTemplate     image.???
//                -start    nSeqStart     1
//                -inc      nSeqIncr      1
//                -num      nNumImages    1000
//                -out      sOutFile      sTemplate(nSeqStart) + "dtav"
//                -outsd    sOutSDFile    sTemplate(nSeqStart) + "dtsd"
//                -avg      sAvgFile      sTemplate(nSeqStart) + "av"
//                -sd       sSDFile       sTemplate(nSeqStart) + "sd"
//                -sigma    fSigma        3.0
//                -max      fMax          1000000.0
//                -min      fMin         -1000000.0
//                -maxsd    fMaxSD        1000000.0
//                -minsd    fMinSD        1.0 (cannot be 0)
//                -noskip                 -skip
//                -warnlimit fWarnlimit   0.001
//                -help
//
//+ToDo
//
//   Error messages need implementing
//
#include "DTMainAverage.h"

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args

{
    CDTMainAverage            oMainAverage;
    
    return  oMainAverage.nExecute(argc, argv);
}
 
