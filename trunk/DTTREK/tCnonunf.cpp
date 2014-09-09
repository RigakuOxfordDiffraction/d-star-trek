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
// tCnonunf.cc        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file tests the classes Cnonunf.
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
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include  "Cstring.h"
#include "Cnonunf.h"

//+Code begin

int main()
{
  Cnonunf     *poNonunf1;
  Cstring      sFiles[3];
  Cimage      *poImage;
  Cstring      sPrefix;

  // Zeroeth check, initialization from a image header

  cout << "Enter name of image file to get nonunf info from: ";
  cin >> sFiles[0];
  poImage = new Cimage(sFiles[0]);

  poImage->oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &sPrefix);
  poImage->oHeader.nGetValue(sPrefix + "NONUNF_INFO", 3, sFiles);
  cout << "3 NONUNF_INFO files are:\n";
  cout << ">>" << sFiles[0] << "<<\n";
  cout << ">>" << sFiles[1] << "<<\n";
  cout << ">>" << sFiles[2] << "<<\n";

  poNonunf1 = new Cnonunf(poImage->oHeader);

  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;


  // First check, initialization to Unknown type

  poNonunf1 = new Cnonunf();
  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;

  // Second check, initialization to Simple type

  cout << "Enter name of simple nonunf image file: ";
  cin >> sFiles[0];
  poNonunf1 = new Cnonunf(sFiles[0]);
  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;

  // Third check, initialization to Dark only type

  cout << "Enter names of 2 nonunf image files: ";
  cin >> sFiles[0] >> sFiles[1];
  poNonunf1 = new Cnonunf(sFiles[0], sFiles[1]);
  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;

  // Fourth check, initialization to Dark etc type

  cout << "Enter names of 3 nonunf image files: ";
  cin >> sFiles[0] >> sFiles[1] >> sFiles[2];
  poNonunf1 = new Cnonunf(sFiles[0], sFiles[1], sFiles[2]);
  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;

  // Fifth check, initialization to Dark etc type, but perhaps errors in images

  cout << "Enter names of 3 nonunf image files: ";
  cin >> sFiles[0] >> sFiles[1] >> sFiles[2];
  poNonunf1 = new Cnonunf(sFiles[0], sFiles[1], sFiles[2]);
  (void) poNonunf1->nList();
  if (poNonunf1->bPixelIsBad(1,1))
    cout << "Pixel 1 1 is bad.\n";
  else
    cout << "Pixel 1 1 is good.\n";
  cout << "Correction factor for 1 1 is " << 
    (poNonunf1->*poNonunf1->m_prfNonunfFactor)(1, 1) << endl;
  cout << "\n\n";
  delete poNonunf1;

  return (0);

}
