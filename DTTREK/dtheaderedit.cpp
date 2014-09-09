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
// dtheaderedit.cc     Initial author: J.W. Pflugrath           18-Apr-1995
//    This edits the header of an image
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
//    dtheaderedit reads from stdin the filename of the input image 
//    to edit the header of and the name of the output image to create.
//    It lists the header to stdout.
//    Next it reads in a keyword value pair  which it adds to the header.
//    Enclose a value with whitespace in single quotes: 'this has whitespace'.
//    The single quotes are removed before adding to the header.
//    When there is an end-of-file on stdin, dtheaderedit writes the output
//    image and exits.
//    dtheaderedit does not change the image in any way, so changing header
//    keywords: DIM, SIZE1, SIZE2, BYTE_ORDER, Data_type, TYPE and
//    COMPRESSION is not allowed.
//
//    Example:  % dtheaderedit 
//                  input.img output.img
//                  NONUNF_FLAG1 1
//                  Description This is the description;
//                  ^D
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "DTMainHeaderEdit.h"

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
    CDTMainHeaderEdit        oMainHeaderEdit;
    
    return  oMainHeaderEdit.nExecute(argc, argv);
}

