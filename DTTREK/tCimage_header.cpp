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
// tCimage_header.cc     Initial author: J.W. Pflugrath           03-Mar-1995
//    This file tests the class Cimage_header.
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
//    This routine tests the class Cimage_header.
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include "Cimage_header.h"
#include "Cstring.h"

int main()
{
  Cimage_header head1;
  Cimage_header head2;
  Cstring sTemp;
  Cstring sEnd = "<<<<\n";
  int ierr;
  int iarr[5] = {0, 1, 2, 3, 5};
  float farr[5] = {0.1, 1.2, 2.3, 3.33, 5.5678901};
  float ffarr[6];
  float f;
  double d;
  int dim1 = 5;
//
// nAddValue functions
//
  cout << "head1 is: " << head1.sGet() << sEnd;

  ierr =  head1.nAddValue("INT", dim1);
  cout << ierr << ": New head1:\n" << head1.sGet() << sEnd;

  ierr = head1.nAddValue("INT", 5);
  cout << ierr << ": Should have error since INT exists\n";
  cout << "head1:\n" << head1.sGet() << sEnd;

  ierr = head1.nAddValue("FLOAT1", (float)3.14157);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("STRING", "this is a string");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nGetValue("STRING", &sTemp);
  cout << ierr << " STRING is>>>" << sTemp << sEnd; 

  ierr = head1.nAddValue("NULL", "");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("FLOAT2", (float) 3.12312313, 2);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("INTWIDE", 3, 5);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("FLOAT3", 3.123123123, 3);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("DOUBLE", 234.24234, 3);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nDelete("FLOAT3");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nDelete("INTWIDE");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nReplaceValue("NULL", 345, 10);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("FLOAT3", 5, farr);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nReplaceValue("FLOAT3", 5, iarr);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nDelete("FLOAT3");
  ierr = head1.nAddValue("FLOAT3", 5, iarr, 10);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nReplaceValue("FLOAT3", 5, farr, 6);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nGetValue("STRING", &sTemp);
  cout << ierr << " STRING is:\n" << sTemp << sEnd; 

  ierr = head1.nGetValue("INT", &dim1);
  cout << ierr << " INT is:\n" << dim1 << sEnd; 

  ierr = head1.nGetValue("FLOAT1", &f);
  cout << ierr << " FLOAT1 is:\n" << f << sEnd; 

  ierr = head1.nGetValue("DOUBLE", &d);
  cout << ierr << " DOUBLE is:\n" << d << sEnd; 

  ierr = head1.nGetValue("FLOAT3", 5, ffarr);
  cout << ierr << " FLOAT3[4] is:\n" << ffarr[4] << sEnd; 

  ierr = head1.nGetValue("FLOAT3", 6, ffarr);
  cout << ierr << " should be error because only 5 elements" << sEnd; 

  sTemp = replicate("abcdefghi  ", 10);
  ierr = head1.nReplaceValue("HEADER_BYTES", sTemp);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  sTemp = "Replace INT with string";
  ierr = head1.nReplaceValue("INT", sTemp);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  sTemp = replicate("123456789 ", 100);
  ierr = head1.nReplaceValue("STRING", sTemp);
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("NOTTEMP1", "This is NOT temp1");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nAddValue("TEMP1", "This is temp1");
  cout << ierr << " head1 is:\n" << head1.sGet() << sEnd; 

  ierr = head1.nGetValue("TEMP1", &sTemp);
  cout << ierr << " TEMP1 is:\n" << sTemp << sEnd; 

  return 0;
}

