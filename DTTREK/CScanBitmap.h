
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CScanBitmap.h            Initial author: Thaddeus Niemeyer      1-17-2002
//  This file contains the member functions of class Cfind which implements
//    the spot finding encapsulation of d*TREK.
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
#include "Cscan.h"
#include "Cimage.h"
#include "dtarray.h"

#ifndef CSCAN_BITMAP_HEADER
#define CSCAN_BITMAP_HEADER

#define CSCAN_BITMAP_HEADERSTRING "dtintegrate bitmap file.\032\012\014"
#define CSCAN_BITMAP_FILE "dtintegrate.dat"

//+Definitions and constants

class DTREK_EXPORT CScanBitmap {

    int m_a2nImageDim[2];           // Image dimensions.
    int m_a2nScanRange[2];          // Images in scan.
    int m_nDataBits;                // Number of data bits for an atypical value.
    int m_nDataBitMask;             // Mask for data bits.
    int m_nTypicalValue;            // Alternate between atypical and typical values in encoding.
    char m_pcHeaderString[50];      // Written to output as a header string.
    itr<int> m_anBitmapStart;       // # of starting bytes from begining.
    Cstring m_sFileName;            // File name.


    
    FILE*   m_pFBitmapDataIn;       // FILE* for input.
    FILE*   m_pFBitmapDataOut;      // FILE* for output.
    int     m_nOutImagePointer;     // Current image that will be written out next. At start equals m_a2nScanRange[0]
    int     m_nInImagePointer;      // Current image that will be read next.  At start equals m_a2nScanRange[0]

    int            m_nLocalBitmapSize;
    unsigned char* m_pcLocalBitmap;
    itr<char>      m_acIOBuffer;            // Buffer to write data into for each image.
    itr<int>       m_anRealToCompressed;    // Maps constructed to map compressed to real and real to compressed values.
    itr<int>       m_anCompressedToReal;


    void  vInit();
    int   nConsistencyCheck();
    int   nWriteHeaderInfo();
    int   nWriteTrailerInfo();
    int   nReadHeaderInfo();
    int   nReadTrailerInfo();

public:

    // basic reading and writing functions.

    int    nOpenOutput(Cstring& sOutput,int nDim0,int nDim1,int nFirstInScan,int nLastInScan,int nDataBits,int nTypicalValue);
    int    nCloseOutput();
    int    nOpenInput(Cstring& sInput);
    int    nOpenInput(Cimage_header& oHeader);  // Only finds the associated bitmap and opens.
    int    nCloseInput();
    int    nWriteNextBitmap();
    int    nReadNextBitmap();
    int    nReadBitmapForScan(Cscan& oScan);    // Gets the bitmap associated with the current pointer for this scan.
    int    nReadBitmapForScan(int nScanNumber);
    
    // read/write the local bitmap object.

    void   vClearLocalBitmap(unsigned char* pcBitmapOther = NULL);                 // Sets to m_nTypicalValue
    void   vSetLocalBitmap(int nPix0,int nPix1,int nValue,unsigned char* pcBitmapOther = NULL);
    void   vSetLocalBitmap(C3Ddata& oShoebox,int nSlice = 0);
    int    nGetLocalBitmap(int nPix0,int nPix1,unsigned char* pcBitmapOther = NULL);
    int    nWriteMaskedImage(Cstring& sImageName,Cimage* poOtherImage);
    int    nMaskOut(Cimage& oImage);
        
    // Update a header to reflect location information for this bitmap object.
    int    nUpdateHeader(Cimage_header& oHeader);

    CScanBitmap();
    ~CScanBitmap();
};

#endif
