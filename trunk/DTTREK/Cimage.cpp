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
// Cimage.cc            Initial author: J.W. Pflugrath           03-Mar-1995
//    This file contains the member functions of class Cimage.
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
//
//   Error messages need implementing
//
//+Include files

#include <string.h>
#include "minmax.h"
#include "dtreksys.h"
#include "Cimage.h"              // Class definition and prototypes
                                 //  Cimage.h includes Cimage_header.h
#include "BitmapRLE.h"
#ifdef WIN32
#include "memory.h"
#endif

#include "dskio.h"
#include "winbmp.h"
#ifdef RAXIS
#include "raxis.h"
#endif
#ifdef SIEMENS
#include "bs.h"
#endif
#ifdef DIP2030
#include "dip2030.h"
#endif
#ifdef MARIP
#include "marip.h"
#endif
#ifdef BAS2000
#include "bas2000.h"
#endif
#ifdef HIPIC
#include "hipic.h"
#endif

#ifdef DTREK_CBF
#include "cbfdtrek.h"
#endif

#if !defined(VC6)
using std::cin;
using std::setw;
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

int Cimage::ms_nVerbose = 1;
int Cimage::ms_nPilatusPixelValueOffset = 0;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cimage::Cimage()
{
  (void) nInitValues();
}

Cimage::Cimage (const int nDim0, const int nDim1)
{

// Construct an image of a given size, but do not allocate space for The_Data

  (void) nInitValues();
  if ( (0 > nDim0) || (0 > nDim1) )
    {
      cout << "     Error constructing image, invalid dimensions!" << endl;
    }
  else
    {
      m_nDim[0]    = nDim0;
      m_nDim[1]    = nDim1;
      m_nData_size = nDataSizeNeeded();

      if (0 < m_nData_size) m_The_Data.pc = new char [m_nData_size];
      m_nBytesAllocated = m_nData_size;
      m_Next_Pixel.pc = m_The_Data.pc;

      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder, sGetByteOrderCPU());
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression,
                                   "None");
      m_eThe_State       = eImage_available_state;
    }
}
/******************
Cimage::Cimage (const int nDim0, const int nDim1, const eImage_data_types nData_type)
{

// Construct an image of a given size and Data_type AND
//   allocate space for The_Data

  (void) nInitValues();
  if ( (0 > nDim0) || (0 > nDim1) )
    {
      cout << "     Error constructing image, invalid dimensions!" << endl;
    }
  else
    {
      m_nDim[0]    = nDim0;
      m_nDim[1]    = nDim1;
      m_eData_type = nData_type;
      m_nData_size = nDataSizeNeeded();

//    m_The_Data.pc = new (char) [m_nData_size];
    if (0 < m_nData_size) m_The_Data.pc = new char [m_nData_size];
    m_nBytesAllocated = m_nData_size;
    m_Next_Pixel.pc = m_The_Data.pc;

    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder, sGetByteOrderCPU());
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression, "None");
  }
}
*************/ 
///////////////////////////////////////////////////////////////////////////////////////////////////////
Cimage::Cimage (const int nDim0, 
                const int nDim1, 
                const eImage_data_types nData_type, 
                const char* sName)
{
   // Construct an image of a given size and Data_type AND
   //   allocate space for The_Data

    (void)nInitValues();

    ////////////////////
    m_sFilename = sName;
    ////////////////////

    if ( (0 > nDim0) || (0 > nDim1) )
    {
        cout << "     Error constructing image, invalid dimensions!" << endl;
        return;
    }
    
    m_nDim[0]    = nDim0;
    m_nDim[1]    = nDim1;
    m_eData_type = nData_type;
    m_nData_size = nDataSizeNeeded();

    if (0 < m_nData_size) m_The_Data.pc = new char [m_nData_size];
    m_nBytesAllocated = m_nData_size;
    m_Next_Pixel.pc = m_The_Data.pc;

    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder, sGetByteOrderCPU());
    (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression, "None");
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

Cimage::Cimage (const int nDim0, const int nDim1,
                const eImage_data_types eData_type,
                void *pvData)
{

// Construct an image of a given size and Data_type AND
//   allocate space for m_The_Data and stuff Data into m_The_Data

  (void) nInitValues();
  if ( (0 > nDim0) || (0 > nDim1) )
    {
      cout << "     Error constructing image, invalid dimensions!" << endl;
    }
  else
    {
      m_nDim[0]    = nDim0;
      m_nDim[1]    = nDim1;
      m_eData_type = eData_type;
      m_nData_size = nDataSizeNeeded();

      if (0 < m_nData_size) m_The_Data.pc     = new char [m_nData_size];
      m_nBytesAllocated = m_nData_size;
      m_Next_Pixel.pc   = m_The_Data.pc;

//      char *pcTo, *pcFrom;
//      pcTo   = m_The_Data.pc;
//      pcFrom = (char *)pvData;
//      for (int i = 0; i < m_nData_size; i++) *pcTo++ = *pcFrom++;

      if (NULL != pvData)
        memcpy(m_The_Data.pv, pvData, m_nData_size);

      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder, sGetByteOrderCPU());
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression, "None");
    }
}

Cimage::Cimage (const Cstring& sName)
{
  int nStat;
  (void) nInitValues();
  nStat = nRead(sName);  // If this fails, what is returned?
  if (0 != nStat)
    {
      if (NULL != m_The_Data.pc)
        {
          delete [] m_The_Data.pc;
          m_The_Data.pc   = NULL;
          m_nBytesAllocated = 0;
        }
      (void) nInitValues();
    }
}

Cimage::Cimage (Cimage& oImageIn, void *pData)
{
  // Construct an image using a Cimage
  // and a pointer to the data which may be NULL
  // The newly constructed image has the same state as the
  // input image even though the data may not have been replaced.

  (void) nInitValues();

  m_eThe_State       = eImage_unknown_state;
  (void) oImageIn.nGetDimensions(&m_nDim[0], &m_nDim[1]);
  m_eData_type       = oImageIn.m_eData_type;
  m_eCompression     = oImageIn.m_eCompression;
  m_nSequence_number = oImageIn.m_nSequence_number;
  m_sFilename        = oImageIn.m_sFilename;
  m_oHeader          = oImageIn.m_oHeader;
  m_fSatVal          = oImageIn.m_fSatVal;
  m_fMinRawPixOKValue= oImageIn.m_fMinRawPixOKValue;
  m_nPixValOffset    = oImageIn.m_nPixValOffset;
  m_nData_size       = oImageIn.m_nData_size;
  m_nBytesAllocated  = 0;
  m_sDescription     = oImageIn.m_sDescription;
  m_sScan_key        = oImageIn.m_sScan_key;
  m_eImgByteOrder    = eCPU_other_endian;
  int i;
  for (i = 0; i < 10; i++) m_nSelection[i] = oImageIn.m_nSelection[i];
  for (i = 0; i < 5; i++)
    {
      m_fCompressInfo[i] = oImageIn.m_fCompressInfo[i];
    }
//+JWP 2007-10-12
  if (1 < m_fCompressInfo[0]) 
    {
      if (eImage_uI2 == m_eData_type)
        {
          prfGetPixel      = &Cimage::fGetRAXISPixel;
          prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
        }
      else
        {
          cout << "SEVERE ERROR in m_eData_type and m_fCompressInfo!\n" << endl;
        }
    }
//-JWP 2007-10-12
  m_eOrientation     = oImageIn.m_eOrientation;

  if (NULL != m_The_Data.pc)
    {
      delete [] m_The_Data.pc;
      m_The_Data.pc = NULL;
    }
  m_nBytesAllocated = nDataSizeNeeded();
  if (-1 == m_nBytesAllocated) m_nBytesAllocated = 0;
  if (0 < m_nBytesAllocated) m_The_Data.pc     = new char [m_nBytesAllocated];

  if (NULL != pData)
    {
      // Copy the data verbatim over

      char *pcTemp1, *pcTemp2;
      pcTemp1 = m_The_Data.pc;
      pcTemp2 = (char *) pData;
      for (i = 0; i < min(m_nBytesAllocated,oImageIn.m_nBytesAllocated); i++)
        {
          *pcTemp1++ = *pcTemp2++;
        }
    }

  if (0 < oImageIn.m_nBitmapMaskByteCount)
  {
      m_nBitmapMaskByteCount = oImageIn.m_nBitmapMaskByteCount;
      
      if (NULL != m_puiBitmapMask)
          delete [] m_puiBitmapMask;
      m_puiBitmapMask = new unsigned short int [m_nBitmapMaskByteCount/2];
      memcpy(m_puiBitmapMask,oImageIn.m_puiBitmapMask,sizeof(unsigned short int)*(m_nBitmapMaskByteCount/2));
  }

  m_eThe_State = oImageIn.m_eThe_State;
}

Cimage::Cimage(C3Ddata *poShoebox, const int nDir, Crefln *poRefln)
{
  // For debugging purposes, create a Cimage from a shoebox
  // so that it can be viewed with the image display program.
  // Place the shoebox data in the image 3 ways: viewed along
  // *poShoebox   pointer to the shoebox object to create an image from
  // nDir         the shoebox plane(s) to place in the 2D image plane
  // *poRefln     pointer to the optional reflection object whose values
  //                 will be placed in the image header.

  // Do some typical image construction stuff.

  // Compute the number of shoebox layers (tiles) in the two image directions
  // Try to make as square an image as possible

  int   nRows, nCols;      // Shoebox layers will be tiled in created image
  float fTemp;

  fTemp  = (float)sqrt((double)poShoebox->nGetExt(2));
  nCols  = (int) (fTemp + 0.5);
  nRows  = poShoebox->nGetExt(2) / nCols;
  if ( (nRows * nCols) < poShoebox->nGetExt(2))
    {
      if (poShoebox->nGetExt(0) > poShoebox->nGetExt(1))
        {
          nRows++;
        }
      else
        {
          nCols++;
        }
    }

  (void) nInitValues();                             // General image init

  m_nDim[0]         = nCols * poShoebox->nGetExt(0);  // Set fast extent
  m_nDim[1]         = nRows * poShoebox->nGetExt(1);  // Set slow extent

  if (e3Ddata_ushort == poShoebox->m_eType)
    {
      m_eData_type  = eImage_uI2;                     // Keep unsigned short
    }
  else
    {
      m_eData_type  = eImage_realIEEE;                // Float data type

      // Re-assign pointer to get float pixel

      prfGetPixel     = &Cimage::fGetfPixel;
    }

  m_eThe_State      = eImage_available_state;         // Be optimistic

  m_nData_size      = nDataSizeNeeded();              // Compute data size
//  m_The_Data.pc     = new (char) [m_nData_size];        // Allocate it
  if (-1 == m_nData_size) m_nData_size = 0;
  if (0 < m_nData_size)
    m_The_Data.pc     = new char [m_nData_size];        // Allocate it
  m_nBytesAllocated = m_nData_size;                     // Save size allocated
  m_Next_Pixel.pc   = m_The_Data.pc;                    // Adjust pointers

  // Fix up header

  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder,
                               sGetByteOrderCPU());
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression, "None");
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDtdisplayOrientation,
                               "-X+Y");
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sShoeboxSize,
                               3, poShoebox->m_nExt);
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sShoeboxOrig, 3,
                               poShoebox->m_nOrigOffset);
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sShoeboxDir, nDir);

  // If available, put refln info into header

  if (NULL != poRefln)
    {
      (void) poRefln->nUpdateHeader(&m_oHeader);
    }

  // Now that an image has been constructed, place the shoebox data
  // into the image data.

  int   i, j, k;           // Loop counters
  int   nPx0, nPx1;        // Pixel coordinate
  int   nRowOff, nColOff;  // Row and column offsets

  nRowOff = 0;
  nColOff = 0;

  // Loop through pixels in the shoebox and place them in the image

  if (eImage_realIEEE == m_eData_type) 
    {
      float *pfTemp;           // A convenience pointer
      double *pflTemp;
      pfTemp = m_The_Data.pf;
      for (k = 0; k < (m_nDim[0] * m_nDim[1]); k++)
        {
            *pfTemp++ = 65536.0;
        }

      if (poShoebox->m_eType == e3Ddata_double)        
        pflTemp = poShoebox->m_pflData;
      else
        pfTemp = poShoebox->m_pfData;

      for (k = 0; k < poShoebox->nGetExt(2); k++)
        {
          // For each layer compute row and column offset
          nPx1 = poShoebox->nGetExt(1) * nRowOff;
          for (j = 0; j < poShoebox->nGetExt(1); j++)
            {
              nPx0 = poShoebox->nGetExt(0) * nColOff;
              for (i = 0; i < poShoebox->nGetExt(0); i++)
                {
              if (poShoebox->m_eType == e3Ddata_double)
                nSetPixel(nPx0++,nPx1,(float) *pflTemp++);
              else
                        nSetPixel(nPx0++, nPx1, *pfTemp++);
                } // end i
              nPx1++;
            }  // end j

          nColOff++;
          if (nColOff >= nCols)
            {
              nColOff = 0;
              nRowOff++;
            }
        }  // end k
    }
  else
    {
      unsigned short int *puiTemp;           // A convenience pointer
      puiTemp = m_The_Data.pui;
      for (k = 0; k < (m_nDim[0] * m_nDim[1]); k++)
        {
          *puiTemp++ = 65535;
        }

      puiTemp = poShoebox->m_puiData;

      for (k = 0; k < poShoebox->nGetExt(2); k++)
        {
          // For each layer compute row and column offset
          nPx1 = poShoebox->nGetExt(1) * nRowOff;
          for (j = 0; j < poShoebox->nGetExt(1); j++)
            {
              nPx0 = poShoebox->nGetExt(0) * nColOff;
              for (i = 0; i < poShoebox->nGetExt(0); i++)
                {
                  nSetPixel(nPx0++, nPx1, *puiTemp++);
                } // end i
              nPx1++;
            }  // end j

          nColOff++;
          if (nColOff >= nCols)
            {
              nColOff = 0;
              nRowOff++;
            }
        }  // end k
    }
}

/****************************************************************************
 ****************************************************************************/
Cimage::Cimage (const Cstring& sHeaderString,
                const void *pData)
{
  int nStat;
  int i;
  Cimage_header *p;

  (void) nInitValues();

  p = new Cimage_header(sHeaderString.string());
  m_oHeader = *p;
  delete p;

  m_eThe_State      = eImage_available_state;         // Be optimistic

  nStat = nSetHeader(NULL);

  if (NULL != m_The_Data.pc)
    {
      delete [] m_The_Data.pc;
      m_The_Data.pc = NULL;
    }
  m_nBytesAllocated = nDataSizeNeeded();
  m_nData_size      = nDataSizeNeeded();              // Compute data size

  if (-1 == m_nData_size) m_nData_size = 0;
  if (0 < m_nData_size)
    m_The_Data.pc     = new char [m_nBytesAllocated];
  m_Next_Pixel.pc   = m_The_Data.pc;                    // Adjust pointers

  if (NULL != pData)
    {
      // Copy the data verbatim over

      char *pcTemp1, *pcTemp2;
      pcTemp1 = m_The_Data.pc;
      pcTemp2 = (char *) pData;
      for (i = 0; i < m_nBytesAllocated; i++)
           {
             *pcTemp1++ = *pcTemp2++;
           }
    }
}

/****************************************************************************
 ****************************************************************************/
Cimage &
Cimage::operator=(const Cimage &oImageIn)
{
  int i;

  // Don't call this because it will cause a memory leak if the image
  // currently has data.  nInitValues() will set the pointer to the
  // data to be NULL, but not deallocate the memory.  
  //(void) nInitValues();

  m_eThe_State = oImageIn.m_eThe_State;
  m_nDim[0] = oImageIn.m_nDim[0];
  m_nDim[1] = oImageIn.m_nDim[1];
  m_eData_type       = oImageIn.m_eData_type;
  m_eCompression     = oImageIn.m_eCompression;
  m_nSequence_number = oImageIn.m_nSequence_number;
  m_sFilename        = oImageIn.m_sFilename;
  m_oHeader          = oImageIn.m_oHeader;
  m_fSatVal          = oImageIn.m_fSatVal;
  m_fMinRawPixOKValue= oImageIn.m_fMinRawPixOKValue;
  m_nData_size       = oImageIn.m_nData_size;
  m_nBytesAllocated  = 0;
  m_sDescription     = oImageIn.m_sDescription;
  m_sScan_key        = oImageIn.m_sScan_key;
  m_eImgByteOrder    = oImageIn.m_eImgByteOrder;
  for (i = 0; i < 10; i++) m_nSelection[i] = oImageIn.m_nSelection[i];
  for (i = 0; i < 5; i++)
    {
      m_fCompressInfo[i] = oImageIn.m_fCompressInfo[i];
    }

  m_eOrientation     = oImageIn.m_eOrientation;

  if (NULL != m_The_Data.pc)
    {
      delete [] m_The_Data.pc;
      m_The_Data.pc = NULL;
    }
  m_nBytesAllocated = nDataSizeNeeded();
  if (-1 == m_nBytesAllocated) m_nBytesAllocated = 0;
  if (0 < m_nBytesAllocated)
    {
      m_The_Data.pc     = new char [m_nBytesAllocated];
      memcpy(m_The_Data.pc,oImageIn.m_The_Data.pc,m_nBytesAllocated);
    }

  m_Next_Pixel.pc = m_The_Data.pc + (oImageIn.m_Next_Pixel.pc-oImageIn.m_The_Data.pc);

  m_nOverflowCount = oImageIn.m_nOverflowCount;
  if(0 != m_nOverflowCount){
                m_pnOverflowOffset = new int   [m_nOverflowCount];
                m_pfOverflowValue  = new float [m_nOverflowCount];
                for (i = 0; i < m_nOverflowCount; i++){
                    m_pnOverflowOffset[i] = oImageIn.m_pnOverflowOffset[i];
                    m_pfOverflowValue[i]  = oImageIn.m_pfOverflowValue[i];
          }
  }
  else{
     m_pnOverflowOffset = NULL;
     m_pfOverflowValue = NULL;
  }

  prfGetPixel = oImageIn.prfGetPixel;
  prfCvtUItoFloat = oImageIn.prfCvtUItoFloat;
  prnSetPixel = oImageIn.prnSetPixel;

  m_nDetNum            = oImageIn.m_nDetNum;
  
  ///////////////////////////////////////////////////////////////////////////
  /// RB 9/15 We need a "deep copy", not a "shallow copy"!
  //m_puiBitmapMask      = oImageIn.m_puiBitmapMask;
  if (NULL != m_puiBitmapMask)
      delete [] m_puiBitmapMask;
          
  m_puiBitmapMask = new unsigned short int [oImageIn.m_nBitmapMaskByteCount/2];
  memcpy(reinterpret_cast<void *>(m_puiBitmapMask),
         reinterpret_cast<void *>(oImageIn.m_puiBitmapMask),
         oImageIn.m_nBitmapMaskByteCount);
  ///////////////////////////////////////////////////////////////////////////

  m_nBitmapMaskByteCount = oImageIn.m_nBitmapMaskByteCount;

  if(NULL != m_pcHeader){
     delete [] m_pcHeader;
     m_pcHeader = NULL;
  }
  if(NULL != oImageIn.m_pcHeader){
     m_pcHeader = new char [sizeof(RAXIS_header)];
     memcpy(reinterpret_cast<void *>(m_pcHeader),
            reinterpret_cast<void *>(m_pcHeader),
            sizeof(RAXIS_header));
  }

  m_poGonioMask = oImageIn.m_poGonioMask;
  m_bReadApplyGonioMask = oImageIn.m_bReadApplyGonioMask;
  m_bReadApplyEmbeddedMask = oImageIn.m_bReadApplyEmbeddedMask;
  m_bWriteRAXIS = oImageIn.m_bWriteRAXIS;

  return *this;
}

Cimage::~Cimage()
{
  if (NULL != m_The_Data.pc)
    {
      delete [] m_The_Data.pc;
      m_The_Data.pc     = NULL;
      m_nBytesAllocated = 0;
    }
  if (NULL != m_pnOverflowOffset)
    {
      delete [] m_pnOverflowOffset;
    }
  if (NULL != m_pfOverflowValue)
    {
      delete [] m_pfOverflowValue;
    }
  if (NULL != m_puiBitmapMask)
    {
      delete [] m_puiBitmapMask;
      m_nBitmapMaskByteCount = 0;
    }
  if (NULL != m_pcHeader)
  {
      delete [] m_pcHeader;
      m_pcHeader = NULL;
  }
  if (NULL != m_poGonioMask)
  {
          delete m_poGonioMask;
          m_poGonioMask = NULL;
  }

  m_nOverflowCount = 0;
}

//+Public functions

int Cimage::nGetDimensions(int *pnDim0, int *pnDim1)
{
  *pnDim0 = m_nDim[0];
  *pnDim1 = m_nDim[1];
  return (m_nDim[0] * m_nDim[1]);
}

// Read an image by prompting for a filename, giving a filename or constructing
//       a filename from a template and sequence number;
int Cimage::nRead()
{
  Cstring sName = "";
  return (nRead(sName));
}

int Cimage::nRead(const Cstring& sName)
{
  int     nStat;
  int     nBytes;
  int     nFile;
  Cstring sTemp;
  Cstring sDetName;
  Cstring sUncompress;
  Cstring sFilename;

  if ("" == sName)
    {
      cout << "Enter the image filename to read: ";
      cin  >> sTemp;
    }
  else
    sTemp = sName;

  sTemp = sTransSymbol(sTemp);

  m_sFilename = sTemp;

//  cout << "     Cimage::nRead filename is " << sTemp << "\n";

  m_eThe_State = eImage_unknown_state;

  // If filename ends in .Z, .gz, or .bz2, then 
  // uncompress the file to a temporary file first,
  // then read it in.

  nBytes    = sTemp.length();
  sFilename = sTemp;
  if ( ('Z' == sTemp.GetAt(nBytes-1)) && ('.' == sTemp.GetAt(nBytes-2)) )
    {
#ifdef WIN32

#ifdef SSI_PC
      //sUncompress = sFileBuildName(GetInstallDirectory(),"gzip");
      //sUncompress += " -d -c ";
      sUncompress = "gzip -d -c ";
#else
      sUncompress = "gzip -d -c ";
#endif  // SSI_PC

#else
      sUncompress = "zcat ";
#endif  // WIN32
      sTemp       = sTemp.before(nBytes-2);
    }
  else if (   ('z' == sTemp.GetAt(nBytes-1))
           && ('.' == sTemp.GetAt(nBytes-3)) && ('g' == sTemp.GetAt(nBytes-2)))
    {
#ifdef WIN32

#ifdef SSI_PC
      //sUncompress = sFileBuildName(GetInstallDirectory(),"gzip");
      //sUncompress += " -d -c ";
      sUncompress = "gzip -d -c ";
#else
      sUncompress = "gzip -d -c ";
#endif  // SSI_PC

#else
      sUncompress = "gunzip -c ";
#endif  // WIN32
      sTemp       = sTemp.before(nBytes-3);
    }
  else if (   ('2' == sTemp.GetAt(nBytes-1)) && ('z' == sTemp.GetAt(nBytes-2))
           && ('b' == sTemp.GetAt(nBytes-3)) && ('.' == sTemp.GetAt(nBytes-4)))
    {
#ifdef WIN32

#ifdef SSI_PC
      //sUncompress = sFileBuildName(GetInstallDirectory(),"bzip2");
      //sUncompress += " -d -c ";
      sUncompress = "";  // Not supported on Windows (TM) for now.
#else
      sUncompress = "bunzip2 -k -c ";
#endif  // SSI_PC

#else
      sUncompress = "bunzip2 -k -c ";
#endif  // WIN32
      sTemp       = sTemp.before(nBytes-4);
    }
  else
    {
      sUncompress = "";
    }

  if ("" != sUncompress)
    {
      // File is most likely compressed, so uncompress it to a temporary file.
      // TODO: probably should put temp file in the /tmp directory

      if ("" != sGetEnv("TMPDIR"))
        {
          // If TMPDIR is defined, then use that as the directory

          sTemp = sTransSymbol("$(TMPDIR)") + '/' + sFileGetBasename(sTemp);
        }

      if (bFileExists(sTemp))
        sTemp = sFileGetTempName(".", "IMG");

//      cout << "...uncompressing " << sFilename << " to " << sTemp << "...\n";
#ifdef SSI_PC
      nDoSystemCommand2(sUncompress + "\"" + sFilename + "\"", sTemp);
#else
      nDoSystemCommand(sUncompress + sFilename + " > " + sTemp);
#endif
//      cout << "...uncompress done.\n";
      sFilename = sTemp;
    }

  // Try to open the file

  nBytes = sFilename.length();
  nFile = 1;
  (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);
  if (0 != nStat)
    {
      cout << "     Error opening file " << sFilename
           << " Error is: " << nStat << "\n";
      if ("" != sUncompress)
        nFileDelete(sFilename);
      return (nStat);
    }
  else
    if (0 < ms_nVerbose)
      {
        if (sFilename == sName)
          {
            cout << "\n     File " << sFilename
                 << " successfully opened.\n";
          }
        else
          {
            cout << "\n     File " << sFilename << "\n     (was "
                 << sName << ") successfully opened.\n";
          }
      }

  //  Next read the header;  What if it is not a d*TREK image???

  nStat = m_oHeader.nRead(&nFile, sFilename, m_sFilename);
  if (0 != nStat)
    {
      cout << "     Error reading header of file " << sFilename
           << "!  Error is: "
           << nStat << "\n";
      (void) dskbcr(&nFile, &nBytes);  // nBytes is temp here
      if ("" != sUncompress)
        nFileDelete(sFilename);
      return (nStat);
    }

  nStat = nSetHeader((Cimage_header*)NULL);

//  Done with header so ...

//  Read the data

//+2010-04-6 JWP      

// Preliminary to allocate enough memory in order to unbin the data in place

  int     nUnbinFactor = 1;
  if ("" != sGetEnv("DTREK_UNBIN"))
    {
      //cout << "DTREK_UNBIN is set A\n";
      // We do not wish to unbin an image that has already been unbinned;
      float a2fUnbinnedValues[2];
      nStat = m_oHeader.nGetValue("DTREK_UNBINNED", 2, &a2fUnbinnedValues[0]);
      if (0 == nStat)
        {
          // Keyword found in header, so already unbinned
          nUnbinFactor = 1;
        }
      else
        nUnbinFactor = 4;
    }
  else
    {
      nUnbinFactor = 1;
    }
//-2010-04-6 JWP

  if (m_eCompression == eImage_no_compression)
    {
      nBytes = nDataSizeNeeded();
      if (0 >= nBytes)
        {
          // Special case, no binary data.  This is legitimate.

          cout << "     WARNING in Cimage::nRead, data size is <= 0.\n";
          (void) dskbcr(&nFile, &nStat);
          if ("" != sUncompress)
            nFileDelete(sFilename);
          if ( (0 == m_nDim[0]) && (0 == m_nDim[1]) )
            {
              m_eThe_State = eImage_available_state;
              return (0);
            }
          else
            {
              return (3);
            }
        }
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor))       // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes)
                m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes)
            m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      if (NULL == m_The_Data.pc)
        {
          cout << "     Error in Cimage::nRead allocating memory!\n";
          nStat = 4;
        }
      sDetName = "";
#ifdef DIP2030
      m_oHeader.nGetValue(Cstring (D_K_DetectorNames), &sDetName);
      if ("DIP2030_" == sDetName) 
        {
          //      cout << "Using fCvtDIP2030toFloat\n";
          prfGetPixel      = &Cimage::fGetDIP2030Pixel;
          prfCvtUItoFloat  = &Cimage::fCvtDIP2030toFloat;
          m_fCompressInfo[0] = 5.0f;
          (void) m_oHeader.nGetValue("DIP2030_Compression", &m_fCompressInfo[0]);
          if (!m_oHeader.bKeywordExists("DIP2030_TRAILER"))
            {
              // Do this ONLY IF trailer exists on the image and has not
              // been read (that is, keyword DIP2030_TRAILER does not exist).

              // Unusual trailer instead of header, so close file and
              // re-open for reading

              (void) dskbcr(&nFile, &nStat);
              nBytes = sFilename.length();
              nFile = 1;
              (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);
              nBytes = nDataSizeNeeded();
            }
        }
#endif
#ifdef MEDOPTICS
      m_oHeader.nGetValue(Cstring (D_K_DetectorNames), &sDetName);
      if (   ( (10 + 1024 * 1024 * 2) == lFileGetSize(m_sFilename) )
               || ( (10 +  512 *  512 * 2) == lFileGetSize(m_sFilename) )
               || ( 45785088               == lFileGetSize(m_sFilename) )
          )
        {
          // A headerless MEDOPTICS or SPring-8 file
          // Close file and re-open for reading
            
          (void) dskbcr(&nFile, &nStat);
          nBytes = sFilename.length();
          nFile = 1;
          (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);
          nBytes = nDataSizeNeeded();
        }
#endif
// Standard d*TREK file fall-through

      // Actually read the contiguous binary data
      
      (void) dskbr(&nFile, m_The_Data.pc, &nBytes, &nStat);
      
      //+ 5-Apr-2002
      // If a bitmap mask occurs at the end of the pixel values, 
      // then read it in.
      
      m_nBitmapMaskByteCount = 0;
      m_oHeader.nGetValue(Cstring(D_K_BitmapSize), &m_nBitmapMaskByteCount);
      if (0 < m_nBitmapMaskByteCount)
        {
          if (NULL != m_puiBitmapMask)
            delete [] m_puiBitmapMask;
          
          m_puiBitmapMask = new unsigned short int [m_nBitmapMaskByteCount/2];
          
          (void) dskbr(&nFile, (char *)m_puiBitmapMask, 
                       &m_nBitmapMaskByteCount, &nStat);
          if (0 != nStat)
            {
              cout << "WARNING:    Error in Cimage::nRead reading bitmap mask!\n";
              delete [] m_puiBitmapMask;
              m_nBitmapMaskByteCount = 0;
            }
        }
      //- 5-Apr-2002      

#ifdef DIP2030
      if (  (0 == nStat) && ("DIP2030_" == sDetName)
            && (!m_oHeader.bKeywordExists("DIP2030_TRAILER")) )
        {
          // Actually read the DIP2030 trailer

          nStat = nReadDIP2030Trailer(nFile, &m_oHeader);
          m_fCompressInfo[0] = 5.0f;
          (void) m_oHeader.nGetValue("DIP2030_Compression", &m_fCompressInfo[0]);
        }
#endif

      //+JWP 02-Feb-2008
      // Some additional code to go with reading a CBF-file t
      // that has been written out in d*TREK format.

      if ( (eImage_I4 == m_eData_type) && ("PILT_" == sDetName) )
        {
          m_fSatVal          = 1048512.0 - 1.0;               //  2**20 - 2**5
          const int nStatus = m_oHeader.nGetValue(Cimage_header::ms_sMinRawPixOKValue, &m_fMinRawPixOKValue);
          if (0 != nStatus) m_fMinRawPixOKValue = float(0);
          m_fCompressInfo[0] = 1.0;
          m_fCompressInfo[1] = 0.0;
          m_nPixValOffset    = ms_nPilatusPixelValueOffset;

          prfGetPixel        = &Cimage::fGetCBFPixel;
          prfCvtUItoFloat    = &Cimage::fCvtCBFtoFloat;
          prnSetPixel        = &Cimage::nSetCBFPixel;
        }
      //-JWP 02-Feb-2008

      if (0 == nStat)
        {
          // m_eImgByteOrder is set in nSetHeader()

          if (m_eImgByteOrder != nGetByteOrderCPU())
            {
              (void) nSwapByteOrder();
              m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
            }
          m_eThe_State = eImage_available_state;
          nStat = 0;               // All successful here
        }
    }  // end of IF-no_compression block
  else if (m_eCompression == eImage_MSC1_compression)
    {
#ifdef RAXIS
      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor))       // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else 
        {
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;

      // Save the binary RAXIS image header if it's available
      if (NULL != m_pcHeader)
        {
          delete [] m_pcHeader;
          m_pcHeader = NULL;
        }
      m_pcHeader = new char[sizeof(RAXIS_header)];
      memcpy(m_pcHeader,&g_oRaxisHeader,sizeof(RAXIS_header));

      if (0 < nBytes)
        nStat = nReadRAXISdata(&nFile, m_The_Data.pui, &m_oHeader);
      if (0 == nStat)
        {
          m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);
          // Set the function pointer to the correct function
          nStat = m_oHeader.nGetValue(Cstring(Cimage_header::ms_sRaxisCompressionRatio),
                                    &m_fCompressInfo[0]);
          if (0 == nStat)
            {
              if ((float)1.0 >= m_fCompressInfo[0])
                {
                  m_fCompressInfo[1] = 0.0;
                }
              else
                {
                  m_fCompressInfo[1] = 32768.0;
                }

              // Look for saturated value in the image header
              
              nStat = m_oHeader.nGetValue(D_K_SaturatedValue, &m_fSatVal);
              if (0 != nStat)
                {
                  // Saturated value is not in the image header so force a value
                  if ((float)8.0 == m_fCompressInfo[0])
                    // It's an raxis2
                    m_fSatVal = 262136.0;                 // 2**18 - 2**3
                  else if ((float)32.0 == m_fCompressInfo[0])
                    // It's an raxis4
                    //+jwp 21-Nov-2002
                    //m_fSatVal = 1048544.0;               //  2**20 - 2**5
                    // Per Katsu Sasaki
                    m_fSatVal = 1048512.0;               //  2**20 - 2**5 - 2**5
                  //-jwp 21-Nov-2002
                  else if ((float)128.0 == m_fCompressInfo[0])
                    m_fSatVal = 4194176.0;
                }
              prfGetPixel      = &Cimage::fGetRAXISPixel;
              prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
            }
          if (m_eImgByteOrder != nGetByteOrderCPU())
            {
              (void) nSwapByteOrder();
              m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
            }
          m_eThe_State = eImage_available_state;
          nStat = 0;               // All successful here
        }
#else
      cout << "     Error in Cimage::nRead: cannot uncompress MSC1 image yet.\n";
      nStat = 5;
#endif
    } // end of IF-MSC1_compression block
  else if (eImage_PCK1_compression == m_eCompression)
    {
#ifdef MARIP
      // MAR345 or otherwise pck'd file

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor))       // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      // Re-open file
      // Read the data

      if (0 < nBytes)
        {
          FILE *fp;

          if ( ( fp = fopen( sFilename.string(), "rb" ) ) == NULL )
            {
              cout << "     Error opening file " << sFilename << endl;
              nStat = 1;
              return (nStat);
            }
//      cout << "About to call get_pck\n" << flush;

          (void) get_pck(fp, m_The_Data.pui);  // How do we know if this worked?
          // Close the file
          fclose (fp);
        }
      nStat = 0;
      if (0 == nStat)
        {
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);

          m_fSatVal = 262136.0;                 // 2**18 - 2**3
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSaturatedValue,
                                         m_fSatVal);
          m_eData_type  = eImage_uI2;       // Default data_type is unsigned short int

          // Set the function pointer to the correct function

          prfGetPixel        = &Cimage::fGetSiemensPixel;
          prfCvtUItoFloat    = &Cimage::fCvtSiemenstoFloat;

          // Fill any overflow information

          m_oHeader.nGetValue("OVERFLOW_COUNT", &m_nOverflowCount);
          if (0 < m_nOverflowCount)
            {
              if (NULL != m_pnOverflowOffset)
                {
                  delete [] m_pnOverflowOffset;
                }
              if (NULL != m_pfOverflowValue)
                {
                  delete [] m_pfOverflowValue;
                }
              m_pnOverflowOffset = new int   [m_nOverflowCount];
              m_pfOverflowValue  = new float [m_nOverflowCount];
              register int i;
              for (i = 0; i < m_nOverflowCount; i++)
                {
                  m_pnOverflowOffset[i] = 0;
                  m_pfOverflowValue[i]  = 0.0;
                }
              (void) m_oHeader.nGetValue("OVERFLOW_VALUES", m_nOverflowCount,
                                         m_pfOverflowValue);
              (void) m_oHeader.nGetValue("OVERFLOW_OFFSET", m_nOverflowCount,
                                         m_pnOverflowOffset);
            }
        }
      if (m_eImgByteOrder != nGetByteOrderCPU())
        {
          (void) nSwapByteOrder();
          m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
        }

      m_eThe_State = eImage_available_state;
      nStat = 0;               // All successful here

      // Convert any pixels values above 32767 into RAXISII style
      // We do this so when we integrate we don't have to mess with the
      // overflow table which is attached to an image.

      m_fCompressInfo[0] = 8.0;
      m_fCompressInfo[1] = 32768.0;

      unsigned short int *puiTemp, iC0, iC1;
      register int i;
      register int k;
      float fTemp;

      iC0 = (short int) m_fCompressInfo[0];
      iC1 = (short int) m_fCompressInfo[1];
      puiTemp = m_The_Data.pui;
      for (i = 0; i < m_nDim[0] * m_nDim[1]; i++, puiTemp++)
        {
          if (32767 < *puiTemp)
            {
              // Value needs R-AXIS style in-place compression

//            cout << "Before: " << *puiTemp;
              *puiTemp = (*puiTemp / iC0) + iC1;
//            cout << " after: " << *puiTemp << " float: "
//                 << fCvtRAXIStoFloat(*puiTemp)
//                   << endl;
            }
        }

      // Replace overflowed values

      for (k = 0; k < m_nOverflowCount; k++)
        {
          if (   (0 <= m_pnOverflowOffset[k])
              && ( (m_nDim[0] * m_nDim[1])  > m_pnOverflowOffset[k])
              && (999999.0 >= m_pfOverflowValue[k])
              && (     0.0 <= m_pfOverflowValue[k]) )
            {
              fTemp = m_pfOverflowValue[k];
              if (m_fSatVal < fTemp) fTemp = m_fSatVal;
              puiTemp = m_The_Data.pui + m_pnOverflowOffset[k];
//            cout << "Before: " << *puiTemp << " overf: " << fTemp << endl;
//            cout << "BeforeB: " << *(puiTemp-1) << " overf: " << fTemp << endl;
//            cout << "BeforeA: " << *(puiTemp+1) << " overf: " << fTemp << endl;

              *puiTemp = (unsigned short int)
                (fTemp / m_fCompressInfo[0] + m_fCompressInfo[1]);
            }
        }

      prfGetPixel      = &Cimage::fGetRAXISPixel;
      prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
#else

      cout << "     Error in Cimage::nRead: cannot uncompress PCK image yet.\n";
      nStat = 5;
#endif
    } // end of IF-PCK1_compression block
  else if (eImage_BS1_compression == m_eCompression)
    {
#ifdef SIEMENS
//+ 21-Apr-2002
// Changes for the new Bruker format as well as better support of the previous format
      m_eData_type  = eImage_uI2;       // Default data_type is unsigned short int
      m_fSatVal = 300000.0;

      // Set the function pointers to the correct function
      // (These may be reset below.)

      prfGetPixel        = &Cimage::fGetSiemensPixel;
      prfCvtUItoFloat    = &Cimage::fCvtSiemenstoFloat;

      // Next line updates m_fSatVal with any environment overrides

      m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &m_fSatVal);

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor) )       // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      // nType will let us know what type of data is found in the image

      int nBSType = 0;

      if (0 < nBytes)
        nStat = nReadSiemensData(&nFile, m_The_Data.pui, &m_oHeader, &nBSType);

      if (0 == nStat)
        {
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);

        }

/************
// Overflow info is now (21-Apr-2002) done inside nReadSiemensData(...)

          // Fill any overflow information

          m_oHeader.nGetValue("OVERFLOW_COUNT", &m_nOverflowCount);
          if (0 < m_nOverflowCount)
            {
                if (NULL != m_pnOverflowOffset)
                  {
                    delete [] m_pnOverflowOffset;
                  }
                if (NULL != m_pfOverflowValue)
                  {
                    delete [] m_pfOverflowValue;
                  }
                m_pnOverflowOffset = new int   [m_nOverflowCount];
                m_pfOverflowValue  = new float [m_nOverflowCount];
                register int i;
                for (i = 0; i < m_nOverflowCount; i++)
                  {
                    m_pnOverflowOffset[i] = 0;
                    m_pfOverflowValue[i]  = 0.0;
                  }
                (void) m_oHeader.nGetValue("OVERFLOW_VALUES", m_nOverflowCount,
                                           m_pfOverflowValue);
                (void) m_oHeader.nGetValue("OVERFLOW_OFFSET", m_nOverflowCount,
                                           m_pnOverflowOffset);
            }
        }

      if (m_eImgByteOrder != nGetByteOrderCPU())
        {
          (void) nSwapByteOrder();
          m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
        }

      // Place overload values into image if value will fit in unsigned short

      register int i;
      for (i = 0; i < m_nOverflowCount; i++)
        {
          if (65535.0 > m_pfOverflowValue[i])
            {
              *(m_The_Data.pui + m_pnOverflowOffset[i])
                = (unsigned short int) m_pfOverflowValue[i];
              m_pfOverflowValue[i] = 0.0;
            }
        }

      // Remove used up overflows from the table

      int nOverflowRemoved = 0;
      for (i = 0; i < m_nOverflowCount; i++)
        {
          if (0.0 == m_pfOverflowValue[i])
            {
              register int j;
              for (j = i; j < m_nOverflowCount-1; j++)
                {
                  m_pfOverflowValue[j]  = m_pfOverflowValue[j+1];
                  m_pnOverflowOffset[j] = m_pnOverflowOffset[j+1];
                }
              nOverflowRemoved++;
            }
        }
      m_nOverflowCount -= nOverflowRemoved;
      m_eThe_State = eImage_available_state;

      nStat = 0;               // All successful here

      if (65535.0 == m_fSatVal)
        {
          // Siemens/Bruker image has max value of 65535, so
          // no need to use special compression functions

          prfGetPixel        = &Cimage::fGetPixel;
          prfCvtUItoFloat    = &Cimage::fCvtUItoFloat;
          prnSetPixel        = &Cimage::nSetPixel;
        }
      else
        {
          // Correlated sum probably used, so values above 65535 are legit
          // Convert any pixels values above 32767 into RAXISII style
          // We do this so when we integrate we don't have to mess with the
          // overflow table which is attached to an image.

          m_fCompressInfo[0] = 8.0;
          m_fCompressInfo[1] = 32768.0;

          unsigned short int *puiTemp, iC0, iC1;
          register int k;
          float fTemp;

          iC0 = (short int) m_fCompressInfo[0];
          iC1 = (short int) m_fCompressInfo[1];
          puiTemp = m_The_Data.pui;
          for (i = 0; i < m_nDim[0] * m_nDim[1]; i++, puiTemp++)
            {
              if (65535 == *puiTemp)
                {
                  // Get value from the overflow table, then needs compression
        
                  fTemp = 65535.0;
                  for (k = 0; k < m_nOverflowCount; k++)
                    {
                      if (m_pnOverflowOffset[k] == i)
                        {
                          fTemp = m_pfOverflowValue[k];
                          break;
                        }
                    }
                  *puiTemp = (unsigned short int)
                              (fTemp / m_fCompressInfo[0] + m_fCompressInfo[1]);
//            cout << "k, m_nOverflowCount, fTemp: " << k << ' '
//              << m_nOverflowCount << ' ' << fTemp << endl;
//            cout << " after: " << *puiTemp << " float: "
//                 << fCvtRAXIStoFloat(*puiTemp)
//                   << endl;
                }
              else if (32767 < *puiTemp)
                {
                  // Value needs R-AXIS style in-place compression

//                  cout << "Before: " << *puiTemp;
                  *puiTemp = (*puiTemp / iC0) + iC1;
//                cout << " after: " << *puiTemp << " float: "
//                 << fCvtRAXIStoFloat(*puiTemp)
//                   << endl;
                }
            }
          m_fSatVal = 262136.0;                 // 2**18 - 2**3
          prfGetPixel      = &Cimage::fGetRAXISPixel;
          prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
        }
***********************/
      m_fSatVal = 262136.0;                 // 2**18 - 2**3
      m_fCompressInfo[0] = 8.0;
      m_fCompressInfo[1] = 32768.0;
      prfGetPixel      = &Cimage::fGetRAXISPixel;
      prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
      if (0 == nStat) m_eThe_State = eImage_available_state;
#else
      cout << "     Error in Cimage::nRead: cannot uncompress BS1 image yet.\n";
      nStat = 5;
#endif
    } // end of IF-BS1_compression block
#ifdef DTREK_CBF
  else if (eImage_CBF_compression == m_eCompression)
    {
      // Try to read the CBF file format

      m_eData_type  = eImage_uI2; // Default data_type will be unsigned short int
      sTemp = Cstring(D_K_UnsignedShortInt);
      (void) m_oHeader.nGetValue(Cstring (D_K_DataType), &sTemp);
      if (Cstring(D_K_LongInt) == sTemp)
        m_eData_type  = eImage_I4; // data_type will be signed long int
      m_fSatVal = 65535.0;

      // Next line updates m_fSatVal with any environment overrides

      m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &m_fSatVal);

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) 
        nBytes = 0;

      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor) )   // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else 
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      //cout << "About to call nReadCBFData\n" << flush;

      if (0 < nBytes)
        nStat = nReadCBFData(sFilename, m_The_Data.pui, &m_oHeader);

      // OK, it's been uncompressed and is no longer in CBF format
      // But it might now be internally in R-AXIS format!

      if (0 == nStat)
        {
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);
        }
      
      //m_eCompression     = eImage_MSC1_compression;

      m_eCompression     = eImage_no_compression;

      if (eImage_uI2 == m_eData_type)
        {
          m_fSatVal          = 1048512.0;               //  2**20 - 2**5
          m_fCompressInfo[0] = 32.0;
          m_fCompressInfo[1] = 32768.0;
          (void) m_oHeader.nReplaceValue(Cstring (D_K_Compression), Cstring(D_K_Raxis));
          (void) m_oHeader.nReplaceValue(Cstring(Cimage_header::ms_sRaxisCompressionRatio),
                                         m_fCompressInfo[0]);

          prnSetPixel      = &Cimage::nSetPixel;
          prfGetPixel      = &Cimage::fGetRAXISPixel;
          prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
        }
      else
        {
          m_fSatVal          = 1048512.0 - 1.0;               //  2**20 - 2**5
          m_fCompressInfo[0] = 1.0;
          m_fCompressInfo[1] = 0.0;
          m_nPixValOffset    = ms_nPilatusPixelValueOffset;

          prfGetPixel        = &Cimage::fGetCBFPixel;
          prfCvtUItoFloat    = &Cimage::fCvtCBFtoFloat;
          prnSetPixel        = &Cimage::nSetCBFPixel;

        }
      if (0 == nStat) m_eThe_State = eImage_available_state;

    } // end of CBF_compression block
#endif
#ifdef BAS2000
  else if (eImage_BAS2000_compression == m_eCompression)
    {
      // Close file and re-read

      (void) dskbcr(&nFile, &nStat);
      nBytes = sFilename.length();
      nFile = 1;
      (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);

      // Try to read the BAS2000 file format

      m_eData_type  = eImage_uI2; // Default data_type will be unsigned short int
      m_fSatVal = 65535.0;

      // Next line updates m_fSatVal with any environment overrides

      m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &m_fSatVal);

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor) )   // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      //cout << "About to call nReadBAS2000Data\n" << flush;
      if (0 < nBytes)
        nStat = nReadBAS2000Data(&nFile, m_The_Data.pui, &m_oHeader);

      // OK, it's been uncompressed and is no longer in BAS2000 format
      // But it is now internally in R-AXIS format!

      if (0 == nStat)
        {
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);
        }
      
      //m_eCompression     = eImage_MSC1_compression;

      m_eCompression     = eImage_no_compression;
      m_fSatVal          = 65535.0;               //  2**20 - 2**5
          int nCCD = 0;
          (void) m_oHeader.nGetValue(Cstring("BAS2000_INF_CCD"), &nCCD);
          if(nCCD == 1) 
              m_fCompressInfo[0] = 0.0;
          else
              m_fCompressInfo[0] = 8.0;

      m_fCompressInfo[1] = 32768.0;
      (void) m_oHeader.nReplaceValue(Cstring (D_K_Compression), Cstring(D_K_Raxis));
      (void) m_oHeader.nReplaceValue(Cstring(Cimage_header::ms_sRaxisCompressionRatio),
                                     m_fCompressInfo[0]);

      prnSetPixel      = &Cimage::nSetPixel;
      prfGetPixel      = &Cimage::fGetRAXISPixel;
      prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
      if (0 == nStat) m_eThe_State = eImage_available_state;
  } //end of BAS2000_compression block
#endif
#ifdef HIPIC
  else if (eImage_HIPIC_compression == m_eCompression)
    {
      // Close file and re-read

      (void) dskbcr(&nFile, &nStat);
      nBytes = sFilename.length();
      nFile = 1;
      (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);

      // Try to read the HiPic file format

      m_eData_type  = eImage_uI2; // Default data_type will be unsigned short int
      m_fSatVal = 65535.0;

      // Next line updates m_fSatVal with any environment overrides

      m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &m_fSatVal);

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor) )   // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes)
                m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      //cout << "About to call nReadHiPicData\n" << flush;
      
      if (0 < nBytes)
        nStat = nReadHiPicData(&nFile, m_The_Data.pui, &m_oHeader);

      // OK, it's been uncompressed and is no longer in HiPic format
      // But it is now internally in R-AXIS format!

      if (0 == nStat)
        {
          (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, m_sFilename);
        }
      
      //m_eCompression     = eImage_MSC1_compression;

      m_eCompression     = eImage_no_compression;
      m_fSatVal          = 65535.0;               //  2**20 - 2**5
      m_fCompressInfo[0] = 0.0;
      m_fCompressInfo[1] = 32768.0;
      (void) m_oHeader.nReplaceValue(Cstring (D_K_Compression), Cstring(D_K_Raxis));
      (void) m_oHeader.nReplaceValue(Cstring(Cimage_header::ms_sRaxisCompressionRatio),
                                     m_fCompressInfo[0]);

      prnSetPixel      = &Cimage::nSetPixel;
      prfGetPixel      = &Cimage::fGetRAXISPixel;
      prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;
      if (0 == nStat) m_eThe_State = eImage_available_state;
  } //end of HIPIC_compression block
#endif
  else if (m_eCompression == eImage_other_compression)
    {
      // Assume it is WinBMP compression
      // Close file and re-read

      (void) dskbcr(&nFile, &nStat);
      nBytes = sFilename.length();
      nFile = 1;
      (void) dskbor(&nFile, sFilename.string(), &nBytes, &nStat);

      m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
      m_eCompression  = eImage_no_compression;
      m_eData_type    = eImage_uI2;

      nBytes = nDataSizeNeeded();
      if (-1 == nBytes) nBytes = 0;
      if (NULL != m_The_Data.pc)
        {
          if (m_nBytesAllocated != (nBytes * nUnbinFactor) ) // Re-use space unless wrong size
            {
              delete [] m_The_Data.pc;
              m_The_Data.pc = NULL;
              if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor]; // pc & pn are the same via union
            }
        }
      else
        {
          if (0 < nBytes) m_The_Data.pc = new char [nBytes * nUnbinFactor];     // Create new space;
        }
      m_nBytesAllocated = nBytes * nUnbinFactor;
      m_nData_size      = nBytes;

      if (0 < nBytes)
        nStat = nReadWindowsBMPData(&nFile, m_The_Data.pui, &m_oHeader);      
      if (0 == nStat)
      {
          m_eThe_State = eImage_available_state;
      }
      // Anything else?
    }
  else
    {
      cout << "     Error in Cimage::nRead, cannot read compressed data of"
           << " Other_type.\n";
      nStat = 5;
    }
#ifdef MARIP
  // TODO: Check for "HIGH" or overflow pixels in MAR case.  If none,
  // you are all set, otherwise you have to do something.
#endif
  // Close the still open file and return

  (void) dskbcr(&nFile, &nBytes);   // nBytes is used as temp status var here

  if ("" != sUncompress)
    nFileDelete(sFilename);

  if (!nStat) 
    nComputeApplyActiveMask();

  // Next line updates m_fSatVal with any environment overrides

  m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &m_fSatVal);

//+2010-04-6 JWP
  if (4 == nUnbinFactor)
    {
      // Unbin the image data and adjust the image header if needed.
      //cout << "DTREK_UNBIN is set C\n";
      nUnbin();
    }
//-2010-04-6 JWP

  return (nStat);
}

int
Cimage::nComputeApplyActiveMask()
{
  if (m_bReadApplyGonioMask)
    {
      int nNumShadows = 0;
      (void) m_oHeader.nGetValue("CRYSTAL_GONIO_NUM_SHADOWS", &nNumShadows); 

      if (0 < nNumShadows)
        {
          CGonioMask oMask;
          // Use the new method of defining shadows:
          oMask.nInitShadows(this);
          oMask.nMaskShadows(this);
          // Ignore the old methods of defining shadows:
/***********************
          // Hard code the values for the detector plate here.
          // This is a TEMPORARY solution;  
          // we would most prefer to have this info in the embedded bitmap.
          // If that's not available, then some other place to get the hard 
          //  coded values would be nice.

          const int nOmegaPoints = 4*3*7;
          double afOmegaPoints[nOmegaPoints] =                          
          { 0.6, -3.1, -0.785, 1.9, -2.6, -0.785,  1.9, -2.6, 0.785, 0.6, -3.1, 0.785, 
            1.9,   -2.6, -0.785, 2.5, -2.0, -0.785,  2.5, -2.0, 0.785, 1.9, -2.6, 0.785, 
            2.5,   -2.0, -0.785, 3.1, -0.8, -0.785,  3.1, -0.8, 0.785, 2.5, -2.0, 0.785,
            3.1,   -0.8, -0.785, 2.9,  1.2, -0.785,  2.9,  1.2, 0.785, 3.1, -0.8, 0.785, 
            2.9,    1.2, -0.785, 2.4,  2.5, -0.785,  2.4,  2.5, 0.785, 2.9,  1.2, 0.785,
            0.6,   -3.1,  0.685, 0.6, -3.1,  1.185,  1.2, -3.1, 1.185, 1.2, -3.1, 0.685,
            0.4,   -3.1, -0.485, 0.6, -3.1, -0.485,  0.6, -3.1, 0.485, 0.4, -3.1, 0.485 
          };                    
          const int nChiPoints = 4*3*3;
          double afChiPoints[nChiPoints] = 
          {
            2.1, -1.05, -0.7, 2.1, 1.05, -0.7, 2.1,  1.05, 0.7, 2.1, -1.05, 0.7,
            1.8,  -0.6, -0.9, 1.8,  0.6, -0.9, 1.8,   0.6, 0.9, 1.8,  -0.6, 0.9,
            1.8,  -0.2, -1.2, 1.8,  0.8, -1.2, 1.8,   0.8, 1.2, 1.8,  -0.2, 1.2
          };
          double* pfPoints;
          const double fScale = 33.0;
          int nPoints;
          int nType;
          itr<double> afX,afY,afZ;
          CGonioMask oMask;
          int nStat;
          int nx,ny;

          oMask.vSetScale(fScale);
                        
          for (nType = 0; nType < 2; nType ++ )
            {
              if (nType == 0)
                {
                  oMask.vSetOmegaStage();
                  pfPoints = &afOmegaPoints[0];
                  nPoints = nOmegaPoints;
                } 
              else 
                {
                  oMask.vSetChiStage();
                  pfPoints = &afChiPoints[0];
                  nPoints = nChiPoints;
                }
              nStat = 0;
              for (ny = 0; (ny < nPoints) && (!nStat); ny+= 12)
                {
                  -afX;
                  -afY;
                  -afZ;
                  for (nx=0; nx<12; )
                    {
                      afX + pfPoints[nx++ + ny];
                      afY + pfPoints[nx++ + ny];
                      afZ + pfPoints[nx++ + ny];
                    }
                  nStat = oMask.nAddBox(afX,afY,afZ);
                }
            }
          if (!nStat)
            nStat = oMask.nReadGonio(m_oHeader);
          if (!nStat)
            nStat = oMask.nBuildMask();
          if (!nStat)
            nStat = oMask.nApplyMask(*this);
          
          if (nStat)
            {
              printf("WARNING:  Could not compute the goniometer mask for RAPID image!\n");
            }
****************************/
        }       
    }

  // The embedded mask should contain active mask info, as well as goniometer mask.
  if (m_bReadApplyEmbeddedMask)
    {
      // See if this image has a bitmap appended to the end.
                
      if (0 < nGetBitmapMaskByteCount())
        {
          BitmapRLE oBitmap;
          if (oBitmap.Set(m_puiBitmapMask,
                          nGetBitmapMaskByteCount()/2))
            {
              if (oBitmap.Apply(*this))
                {
                }
              else
                cout << "ERROR: Embedded bitmap Mask found, but unusable.\n" << flush;
            }
          else
            {
              // ERROR
              cout << "ERROR: Embedded bitmap Mask found, but unusable.\n" << flush;
            }
        }
    }
  return 0;
}

int Cimage::nWrite(void)
{
  Cstring sName;
  cout << "Enter the image filename to write: ";
  cin  >> sName;
  cout << "     Filename to write: " << sName << "\n";
  return (nWrite(sName));
}

int Cimage::nWrite(const Cstring& sName, bool bMustWriteDTrekHeader)
{
  int     nStat;
  int     nBytes;
  int     nFile;
  int     nSize;

  int n;
  Cstring sTemp;
  Cstring sFilenameWrite;

  // Try to open the file

  if ( (NULL == m_The_Data.pc) && (0 !=m_nDim[0]) && (0 != m_nDim[1]) )
    {
      cout << "     Error in Cimage::nWrite, no data available to write!\n";
      return (1);
    }

  sFilenameWrite = sTransSymbol(sName);

  // If the overwrite environment variable is set, then make sure
  // any existing image file gets a version number appended, so it is
  // not overwritten

  nStat = nFileAppendVersion(sFilenameWrite, TRUE);
  if (0 != nStat)
    {
      // Images are so important, this is a FATAL error

      cout << "FATAL ERROR renaming image file: " << sTemp << "!\n";
#ifdef SSI_PC
          return nStat;
#else
      exit (nStat);
#endif
    }

#ifdef RAXIS
  //+2012-09-12 JWP
  // If special environment variable is set,  then ...
  //   Do NOT Write R-AXIS binary image format, instead write d*TREK format with ASCII header
  if ("" != sGetEnv("DTREK_RAXIS_WRITE_ASCII"))
    // Note how the preceding if statement has the if statement below as the block ...
    //-2012-09-12 JWP

  if ((eImage_MSC1_compression == m_eCompression) && (m_bWriteRAXIS) && !bMustWriteDTrekHeader)
    {
      // Write only if there are pixels in this image
      if (0 < nDataSizeNeeded())
        return (nWriteRAXIS(*this, sFilenameWrite,(_RAXIS_header*) m_pcHeader));
    }
#endif

  // programming logic hypothesis (mrp): write to temporary before renaming to desired sFilenameWrite 
  // so that outside consumers of this file won't find it until it's ready to be read.
  Cstring sTempTemp;
  sTempTemp = sFilenameWrite + "TMP";
  //  cout << "A. sTempTemp>>" << sTempTemp << "<<\n" << endl;
  nStat = nFileAppendVersion(sTempTemp, TRUE);
  //  cout << "B. sTempTemp>>" << sTempTemp << "<<\n" << endl;
  nBytes = sTempTemp.length();
  nFile = 1;
  nSize = nDataSizeNeeded();
  if (-1 == nSize) nSize = 0;
  nSize = m_oHeader.nLength() + nSize;

//  (void) dskbow(&nFile, (const char *)sTemp, &nBytes, &nSize, &nStat);
  (void) dskbow(&nFile, sTempTemp.string(), &nBytes, &nSize, &nStat);

  if (0 != nStat)
    {
      cout << "     Error opening file " << sTempTemp << "Error is: " << nStat << "\n";
    return (nStat);
  }
  else
    if (0 < ms_nVerbose)
      {
        cout << "     Temporary file " << sTempTemp << " successfully opened.\n";
      }

  m_sFilename         = sName;                     // Keep any environment var
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sFilename, sName); // stuff

  //  Next update the header (so KEYWORDS really say what will be in the file!)
  //  Data_type, COMPRESSION, BYTE_ORDER

  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sByteOrder, sGetByteOrderCPU());
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sDataType, sGetDataType());

  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSaturatedValue, m_fSatVal);

  // For backwards compatibility, if keyword "TYPE" exists, it must be
  // in the header AFTER "Data_type" keyword
  if (0 == m_oHeader.nGetValue(Cimage_header::ms_sType, &sTemp))
    {
      (void) m_oHeader.nReplaceValue(Cimage_header::ms_sType, sTemp);
    }

  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sCompression, "None");
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, m_nDim[0]);
  (void) m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, m_nDim[1]);

  //  Next write the header

  nStat = m_oHeader.nWrite(&nFile);
  if (0 != nStat)
    {
      cout << "     Error writing header of file " << sName << "!  Error is: "
           << nStat << "\n";
      (void) dskbcw(&nFile, &n);

          // Arbitrary pause to let the operating system tidy up the file before attempting to rename it
          nWait(100);

      if (0 != nFileRename(sTempTemp, sFilenameWrite))
        {
          cout << "     Error renaming file " << sTempTemp << " to " 
               << sFilenameWrite
               << "!  Error is: "
           << nStat << "\n";
        }
      return (nStat);
    }
  nBytes = nDataSizeNeeded();

  if (0 < nBytes)
    {
      (void) dskbw(&nFile, m_The_Data.pc, &nBytes, &nStat);
      if ( (0 == nStat) && (0 < m_nBitmapMaskByteCount) )
        {
          (void) dskbw(&nFile, (char*)m_puiBitmapMask, &m_nBitmapMaskByteCount, &nStat);
          if (0 != nStat)
            {
              cout << "     Error in Cimage::nWrite writing bitmap mask!\n";
              nStat = 5;
            }
        }

      if (0 != nStat)
        {
          cout << "     Error writing file " << sName << "!  Error is: "
               << nStat << "\n";
        }
      else
        {
          if (0 < ms_nVerbose)
            {
              cout << "     Success writing file " << sName << "!\n";
            }
        }
    }
  (void) dskbcw(&nFile, &nStat);

  // Arbitrary pause to let the operating system tidy up the file before attempting to rename it
  nWait(100);

  if (0 != nFileRename(sTempTemp, sFilenameWrite))
    {
      cout << "     Error renaming file " << sTempTemp << " to " 
           << sFilenameWrite
           << " !\n" << flush;
    }
  return (nStat);
}

int Cimage::nList(const int nRect[4])
{

  //  List a rectangle of data specified by nRect (center and widths)
  // To Do:  Format for other data types and orientations

  int nStat;
  int nFast, nSlow;
  int nTemp[4];

  // Adjust request to fit available data

  nTemp[0] = max(0, nRect[0]-(nRect[2]/2)); // (a>b)?a:b;
  nTemp[1] = max(0, nRect[1]-(nRect[3]/2));
  nTemp[0] = min(m_nDim[0]-1, nTemp[0]);
  nTemp[1] = min(m_nDim[1]-1, nTemp[1]);
  nTemp[2] = min(m_nDim[0]-1, nRect[0]+(nRect[2]/2));
  nTemp[3] = min(m_nDim[1]-1, nRect[1]+(nRect[3]/2));
  nTemp[2] = max(0, nTemp[2]);
  nTemp[3] = max(0, nTemp[3]);
  nTemp[2] = nTemp[2] - nTemp[0] + 1;
  nTemp[3] = nTemp[3] - nTemp[1] + 1;

  short int *paInt;
  if (nTemp[2] > 0 && nTemp[3] > 0)
    // Get enough memory for the Rect

//    paInt = new (short int) [nTemp[2]*nTemp[3]];
    paInt = new short [nTemp[2]*nTemp[3]];
  else {
    cout << "     Error Cimage::nList: bad new [" << nTemp[2]*nTemp[3] <<"]!\n";
    return (1);
  }

  nStat = nGetRect(nTemp, (char **)&paInt);
  if (0 != nStat)
    {
      cout << "     Error Cimage::nList: bad nGetRect,nStat = " << nStat << "!\n";
      return (nStat);
    }

  float fValue;
  short int *pLoop;
  pLoop = paInt;
  cout << "\n   Fast>";
  for (nFast = nTemp[0]; nFast < nTemp[0]+nTemp[2]; nFast++) // fast
      cout << setw(6) << nFast;
  cout << "\nSlow\\/";
  for (nSlow = nTemp[1]; nSlow < nTemp[1]+nTemp[3]; nSlow++)
    {
      // slow
      cout << endl << setw(6) << nSlow << ": ";
      for (nFast = nTemp[0]; nFast < nTemp[0]+nTemp[2]; nFast++)
        { // fast
          fValue = (this->*prfCvtUItoFloat)(*pLoop++);
            cout << setw(6) << (int)fValue;
          //cout << setw(6) << *pLoop++;
        }
      cout << "  :" << nSlow;
    }
  cout << "\nSlow/\\" << endl;
  cout << "   Fast>";
  for (nFast = nTemp[0]; nFast < nTemp[0]+nTemp[2]; nFast++) // fast
      cout << setw(6) << nFast;
  cout << "\n\n";

  delete [] paInt;
  paInt = NULL;

  return (0);
}

float
Cimage::fGetfPixel(const int nWhere0, const int nWhere1)
{
    // Get a pixel at location nWhere[0], nWhere[1]

    // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (0.0);
#endif
  return ( *(m_The_Data.pf+nWhere1*m_nDim[0] + nWhere0));
}

float
Cimage::fGetuiPixel(const int nWhere0, const int nWhere1)
{
    // Get a pixel at location nWhere[0], nWhere[1]

    // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (0.0);
#endif
  return ( (float) *(m_The_Data.pui+nWhere1*m_nDim[0] + nWhere0));
}

float Cimage::fGetPixel(const int nWhere0, const int nWhere1)
{
  // Get a pixel at location nWhere[0], nWhere[1]

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (0.0);
#endif
  if (m_eData_type == eImage_I2)
    return ( (float) *(m_The_Data.pi+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_byte)
    return ( (float) *(m_The_Data.pc+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_ubyte)
    return ( (float) *(m_The_Data.puc+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_uI2)
    return ( (float) *(m_The_Data.pui+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_I4)
    return ( (float) *(m_The_Data.pl+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_uI4)
    return ( (float) *(m_The_Data.pul+nWhere1*m_nDim[0] + nWhere0));
  else if (m_eData_type == eImage_realIEEE)
    return ( *(m_The_Data.pf+nWhere1*m_nDim[0] + nWhere0));
  else
    return (0.0);
}

int Cimage::nGetPixel(const int nWhere[2], short int *piPixel)
{
  // Get a pixel at location nWhere[0], nWhere[1]

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere[0])
      || (0 > nWhere[1])
      || (nWhere[0] >= m_nDim[0])
      || (nWhere[1] >= m_nDim[1]) )
    return (1);
#endif
  *piPixel = *(m_The_Data.pi+nWhere[1]*m_nDim[0] + nWhere[0]);
  return (0);
}

int Cimage::nGetPixel(const int nWhere0, const int nWhere1, short int *piPixel)
{
  // Get a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *piPixel = *(m_The_Data.pi +  nWhere1 * m_nDim[0] + nWhere0);
  return (0);
}

int Cimage::nGetPixel(const int nWhere0, const int nWhere1,
                      unsigned short int *puiPixel)
{
  // Get a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *puiPixel = *(m_The_Data.pui +  nWhere1 * m_nDim[0] + nWhere0);
  return (0);
}

int Cimage::nGetPixel(const int nWhere[2], int *pnPixel)
{
    // Get a pixel at location nWhere[0], nWhere[1]

    // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere[0])
      || (0 > nWhere[1])
      || (nWhere[0] >= m_nDim[0])
      || (nWhere[1] >= m_nDim[1]) )
    return (1);
#endif
  *pnPixel = *(m_The_Data.pn+nWhere[1]*m_nDim[0] + nWhere[0]);
  return (0);
}

int Cimage::nGetPixel(const int nWhere0, const int nWhere1, int *pnPixel)
{
  // Get a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *pnPixel = *(m_The_Data.pn +  nWhere1 * m_nDim[0] + nWhere0);
  return (0);
}


int Cimage::nGetPixel(const int nWhere0, const int nWhere1, LONG *plPixel)
{
    // Get a pixel at location nWhere0, nWhere1

    // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *plPixel = *(m_The_Data.pl +  nWhere1 * m_nDim[0]  + nWhere0);
  return (0);
}

int Cimage::nGetPixel(const int nWhere0, const int nWhere1,
                      unsigned int *punPixel)
{
  // Get a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *punPixel = *(m_The_Data.pun +  nWhere1 * m_nDim[0] + nWhere0);
  return (0);
}

int Cimage::nSetNextPixel(const int nWhere[2])
{
  return (nSetNextPixel(nWhere[0], nWhere[1]));
}

int Cimage::nSetNextPixel(const int nWhere0, const int nWhere1)
{
  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  m_Next_Pixel.pc = m_The_Data.pc + ((nWhere1*m_nDim[0] + nWhere0)
                                 * nBytesPerPixel());
  return (0);
}

int Cimage::nSetPixel(const int nWhere[2], const short int iPixel)
{
  // Set a pixel at location nWhere[0], nWhere[1]

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere[0])
      || (0 > nWhere[1])
      || (nWhere[0] >= m_nDim[0])
      || (nWhere[1] >= m_nDim[1]) )
    return (1);
#endif
  *(m_The_Data.pi + nWhere[1]*m_nDim[0] + nWhere[0]) = iPixel;

  return (0);
}


int Cimage::nSetPixel(const int nWhere0 , const int nWhere1, const LONG lPixel)
{
  // Set a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *(m_The_Data.pl + nWhere1*m_nDim[0] + nWhere0) = lPixel;

  return (0);
}


int Cimage::nSetPixel(const int nWhere0, const int nWhere1,
                      const unsigned short int uiPixel)
{
  // Set a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *(m_The_Data.pui +  nWhere1 * m_nDim[0] + nWhere0) = uiPixel;

  return (0);
}

int Cimage::nSetPixel(const int nWhere0, const int nWhere1,
                      const ULONG ulPixel)
{
  // Set a pixel at location nWhere0, nWhere1

  // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *(m_The_Data.pul +  nWhere1 * m_nDim[0] + nWhere0) = ulPixel;

  return (0);
}

int Cimage::nSetPixel(const int nWhere0, const int nWhere1,
                      const short int iPixel)
{
    // Set a pixel at location nWhere0, nWhere1

    // Bounds checking... this will really slow us down!
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif
  *(m_The_Data.pi +  nWhere1 * m_nDim[0] + nWhere0) = iPixel;
  return (0);
}

int Cimage::nSetPixel(const int nWhere0, const int nWhere1,
                      const float fPixel)
{
  // Set a pixel at location nWhere0, nWhere1
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  // Bounds checking... this will really slow us down!
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif

  if (m_eData_type == eImage_I2)
    *(m_The_Data.pi+nWhere1*m_nDim[0] + nWhere0)  = (short int) fPixel;
  else if (m_eData_type == eImage_uI2)
    {
#ifdef RAXIS
      if (eImage_MSC1_compression != m_eCompression && m_fCompressInfo[0] <= 1.0f )
        {
#endif
          *(m_The_Data.pui+nWhere1*m_nDim[0] + nWhere0)
            = (unsigned short int) fPixel;
#ifdef RAXIS
        }
      else
        {
          // RAXIS style compression

          if (fPixel <= 32767.0)
            {
              *(m_The_Data.pui+nWhere1*m_nDim[0] + nWhere0)
                = (unsigned short int) fPixel;
            }
          else
            {
              *(m_The_Data.pui+nWhere1*m_nDim[0] + nWhere0)
                = (unsigned short int) ((fPixel / m_fCompressInfo[0]) +
                                        m_fCompressInfo[1]);
            }
        }
#endif
    }
  else if (m_eData_type == eImage_byte)
    *(m_The_Data.pc+nWhere1*m_nDim[0] + nWhere0)  = (char) fPixel;
  else if (m_eData_type == eImage_ubyte)
    *(m_The_Data.puc+nWhere1*m_nDim[0] + nWhere0) = (unsigned char) fPixel;
  else if (m_eData_type == eImage_I4)
    *(m_The_Data.pl+nWhere1*m_nDim[0] + nWhere0)  = (LONG) fPixel;
  else if (m_eData_type == eImage_uI4)
    *(m_The_Data.pul+nWhere1*m_nDim[0] + nWhere0) = (ULONG) fPixel;
  else if (m_eData_type == eImage_realIEEE)
    *(m_The_Data.pf+nWhere1*m_nDim[0] + nWhere0)  = fPixel;
  else
    return (2);

  return (0);
}

int Cimage::nSetCBFPixel(const int nWhere0, const int nWhere1,
                      const float fPixel)
{
  // Set a pixel at location nWhere0, nWhere1
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  // Bounds checking... this will really slow us down!
  if (   (0 > nWhere0)
      || (0 > nWhere1)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (1);
#endif

  *(m_The_Data.pn+nWhere1*m_nDim[0] + nWhere0) = (int) fPixel + m_nPixValOffset;
  return (0);
}

int Cimage::nGetRect(const int nRect[4], char **ppaData)
{
    // Allocate memory for, then get a rectangle of data out of an image and
    //   return a pointer to the data.

  int nFast, nSlow;
  int nStart;

  // Do some error checking.

  if (NULL == m_The_Data.pc)
    // No data available
    return (1);

  if (m_eData_type == eImage_compressed || m_eData_type == eImage_other_type)
    // Unimplemented data_type error:
    return (2);

  if ( (nRect[0] < 0) ||
       (nRect[1] < 0) ||
       (nRect[2] < 1) ||
       (nRect[3] < 1) ||
       ( (nRect[0]+nRect[2]) > m_nDim[0]) ||
       ( (nRect[1]+nRect[3]) > m_nDim[1]) )
    // Out of bounds, invalid arguments
    return (3);

  if (NULL != *ppaData)
    {
      // Allocate memory for result, otherwise use location passed to this array.
      int nBytes = nRect[2] * nRect[3] * nBytesPerPixel();
      if (nBytes <= 0)
        // Another way of finding if no data;
        return (5);
      //    *ppaData = new (char) [nBytes];
      *ppaData = new char [nBytes];
    }

  // Extract the rectangle of data from the image and place in ppaData

  if (m_eData_type == eImage_byte)
    return (6);
  else if (m_eData_type == eImage_ubyte)
    return (6);
  else if (m_eData_type == eImage_I2)
    {
      short int *piTo;
      piTo  = (short int *) *ppaData;
      short int *piFrom;

//            slowstart * fastdim   + faststart
      nStart =  nRect[1] * m_nDim[0]  +  nRect[0];
      for (nSlow = 0; nSlow < nRect[3]; nSlow++)
        {
          // Slow
          piFrom = m_The_Data.pi + nStart;
          for (nFast = 0; nFast < nRect[2]; nFast++)      // Fast
            *piTo++ = *piFrom++;         // Extract from input and stuff in output
          nStart = nStart + m_nDim[0];     // m_nDim[0] is fast dimension of array
        }
      return (0);
    }
  else if (m_eData_type == eImage_I4)
    return (6);
  else if (m_eData_type == eImage_uI2)
    {
      unsigned short int *puiTo;
      puiTo  = (unsigned short int *) *ppaData;
      unsigned short int *puiFrom;

//              slowstart * fastdim   + faststart
      nStart =  nRect[1] * m_nDim[0]  +  nRect[0];
      for (nSlow = 0; nSlow < nRect[3]; nSlow++)
        {      // Slow
          puiFrom = m_The_Data.pui + nStart;
          for (nFast = 0; nFast < nRect[2]; nFast++)      // Fast
            *puiTo++ = *puiFrom++;       // Extract from input and stuff in output
          nStart = nStart + m_nDim[0];     // m_nDim[0] is fast dimension of array
        }
      return (0);
    }
  else if (m_eData_type == eImage_uI4)
    return (6);
  else if (m_eData_type == eImage_realIEEE)
    return (6);

  return (0);
}

int Cimage::nCopyRect(const Cimage& hFrom, const int nFrom[4], const int nTo[2])
{

  //  Get a rectangle of data out of an image and
  //   put in the pre-allocated pointer to the data.

  int nFast, nSlow;
  int nStartFrom;
  int nStartTo; 
        
  // Do some error checking.
  // Bounds checking and will it fit?

  /**
  if (   (2 != hFrom.nBytesPerPixel())
      || (2 !=       nBytesPerPixel()) )
    {
      cout << "ERROR! Cimage::nCopyRect() supports only 2-byte datatypes!\n"
           << flush;
      return (2);
    }
  **/
  if (hFrom.m_eData_type != m_eData_type)
    {
      // Mismatched data types, in the future we'll fix this
      // If bytes per pixel is same, then don't worry about it

      if ((const int)hFrom.nBytesPerPixel() != (const int)nBytesPerPixel())
        {
          cout << "ERROR: Mismatched bytes per pixel!\n" << flush;
          return (2);
        }
      else
        {
          cout << "WARNING: Mismatched image data types, IGNORED!\n" << flush;
        }
    }

  if ( (nFrom[0] < 0) ||                        // Fast px origin before image
       (nFrom[1] < 0) ||                        // Slow px origin before image
       (nFrom[2] < 1) ||                        // Fast extent bogus
       (nFrom[3] < 1) ||                        // Slow extent bogus
       (nFrom[0] > hFrom.m_nDim[0]) ||            // Fast px origin past image
       (nFrom[1] > hFrom.m_nDim[1]) ||            // Slow px origin past image
       (nFrom[0]+nFrom[2] > hFrom.m_nDim[0]) ||   // Fast extent won't fit
       (nFrom[1]+nFrom[3] > hFrom.m_nDim[1]) ||   // Slow extent won't fit
       (nTo[0] < 0) ||                          // Fast px To: before image
       (nTo[1] < 0) ||                          // Slow px To: before image
       ( (nTo[0]+nFrom[2]) > m_nDim[0]) ||        // Slow extent To: won't fit
       ( (nTo[1]+nFrom[3]) > m_nDim[1]) )         // Slow extent To: won't fit
    return (3);                        // Out of bounds, invalid arguments

  if (hFrom.m_The_Data.pi == NULL || m_The_Data.pi == NULL)
    {
      // No data to copy from or copy to
      return (5);
    }

  //+2013-01-07 JWP
  if (2 ==  (const int)hFrom.nBytesPerPixel())
    {
      // Extract the rectangle of data from the image and place in ppaData

      short int *piTo;
      short int *piFrom;

      //            slowstart * fastdim   +         faststart

      nStartFrom =  nFrom[1] * hFrom.m_nDim[0]  +  nFrom[0];
      nStartTo   =    nTo[1] *       m_nDim[0]  +    nTo[0];

      for (nSlow = 0; nSlow < nFrom[3]; nSlow++)
        {      // Slow
          piFrom = hFrom.m_The_Data.pi + nStartFrom;
          piTo   =       m_The_Data.pi + nStartTo;
          for (nFast = 0; nFast < nFrom[2]; nFast++)     // Fast
            *piTo++ = *piFrom++;          // Extract from input and stuff in output
          nStartFrom = nStartFrom + hFrom.m_nDim[0];  // m_nDim[0] is fast dim
          nStartTo   = nStartTo   +       m_nDim[0];  // m_nDim[0] is fast dim
        }
    }
  else if (4 ==  (const int)hFrom.nBytesPerPixel())
    {
      // Extract the rectangle of data from the image and place in ppaData

      int *pnTo;
      int *pnFrom;

      //            slowstart * fastdim   +         faststart

      nStartFrom =  nFrom[1] * hFrom.m_nDim[0]  +  nFrom[0];
      nStartTo   =    nTo[1] *       m_nDim[0]  +    nTo[0];

      for (nSlow = 0; nSlow < nFrom[3]; nSlow++)
        {      // Slow
          pnFrom = hFrom.m_The_Data.pn + nStartFrom;
          pnTo   =       m_The_Data.pn + nStartTo;
          for (nFast = 0; nFast < nFrom[2]; nFast++)     // Fast
            *pnTo++ = *pnFrom++;          // Extract from input and stuff in output
          nStartFrom = nStartFrom + hFrom.m_nDim[0];  // m_nDim[0] is fast dim
          nStartTo   = nStartTo   +       m_nDim[0];  // m_nDim[0] is fast dim
        }
    }
  else
    {
      cout << "ERROR: Unsupported bytes per pixel!\n" << flush;
      return (3);
    }
  //-2013-01-07 JWP
  return (0);
}

//
//+Private functions

// Initialize initial values

int Cimage::nInitValues(void)
{

   // If you make a change here, also make it in the operator=()
   // as it does not call this routine.

  m_eThe_State         = eImage_unknown_state;
  m_eImgByteOrder      = eCPU_other_endian;
  m_nDim[0]            = 0;
  m_nDim[1]            = 0;
  m_eData_type         = eImage_uI2;
  m_eCompression       = eImage_no_compression;
  m_nSequence_number   = -1;
  m_sFilename          = "Internally represented";
  m_fSatVal            = 65535;                     // We don't know this yet;
  m_fMinRawPixOKValue  = float(1);
  m_sDescription       = "Just initialized image";

  m_nDetNum            = 0;
  m_nData_size         = 0;
  m_nBytesAllocated    = 0;
  m_nPixValOffset      = 0;
  m_The_Data.pc        = NULL;
  m_Next_Pixel.pc      = NULL;

  m_puiBitmapMask      = NULL;
  m_nBitmapMaskByteCount = 0;

  m_pcHeader           = NULL;

  // Initialize overflow table to have no overflows

  m_nOverflowCount   = 0;
  m_pnOverflowOffset = NULL;
  m_pfOverflowValue  = NULL;

  // Set default pointers to member functions

  prfGetPixel        = &Cimage::fGetPixel;
  prfCvtUItoFloat    = &Cimage::fCvtUItoFloat;
  prnSetPixel        = &Cimage::nSetPixel;

  if (NULL != m_The_Data.pc) cout << "WARNING WARNING WARNING!\n";
  m_sScan_key          = "unknown";
  int i; for (i = 0; i < 10; i++) m_nSelection[i] = 0;

  for (i = 0; i < 5; i++) m_fCompressInfo[i] = 1.0;

  m_eOrientation        = eO_unknown;

  m_poGonioMask = NULL;
  m_bReadApplyGonioMask = false;
  m_bReadApplyEmbeddedMask = false;
  m_bWriteRAXIS = true;
  
  if ("" != sGetEnv("DTREK_IMAGE_APPLYMASK"))
  {
      m_bReadApplyGonioMask = TRUE;
      m_bReadApplyEmbeddedMask = TRUE;
  }
  
  // Add two things here for backwards compatibility
  m_oHeader.nAddValue(Cimage_header::ms_sType, D_K_MadMarty);
  m_oHeader.nAddValue(Cimage_header::ms_sDim, (int) 2);

  m_nRef = 0;

  for(int ii=0; ii < 5; ii++)
      m_fCompressInfo[ii] = 1.0f;

  return 0;
}


int Cimage::nSwapByteOrder(void)
{
  int i;
  int nBytes = nDataSizeNeeded();
  if ( (m_The_Data.pc == NULL) || (nBytes <=0) || (m_nBytesAllocated == 0) )
    {
      cout << "     Error in Cimage::nSwapByteOrder: no data to swap!\n.";
      return -1;
    }

  if (2 == nBytesPerPixel())
    {
      unsigned short int *puiTemp;
      puiTemp = m_The_Data.pui;
//      cout << "     Info in Cimage::nSwapByteOrder: swapping short int.\n";
      //      for (i = 0; i < nBytes/2; i++)
      i = nBytes / 2;
      while (0 < i)  // Tests show this is 10-20% faster than the for loop.
        {
          *puiTemp = ((*puiTemp & 0xff) << 8) | ((*puiTemp) >> 8);
          puiTemp++;
          i--;
        } 
      return 0;
    }
  else if (4 == nBytesPerPixel())
    {
      unsigned int *punTemp;
      punTemp = m_The_Data.pun;
      i = nBytes / 4;
      while (0 < i)
        {
          *punTemp =   (*punTemp << 24) | ((*punTemp << 8) & 0x00ff0000)
                      | ((*punTemp >> 8) & 0x0000ff00) | (*punTemp >> 24);
          punTemp++;
          i--;
        }
      return 0;
    }
  else if (1 == nBytesPerPixel())
    return 0;
  else
    {
      cout << "     Error in CIMAGE::nSwapByteOrder unknown bytes\n.";
      return 2;
    }
}


int Cimage::nBytesPerPixel(void) const
{
  if (m_eData_type == eImage_byte || m_eData_type == eImage_ubyte)
    return sizeof(char);
  else if (m_eData_type == eImage_I2 || m_eData_type == eImage_uI2)
    return sizeof(short int);
  else if (m_eData_type == eImage_I4 || m_eData_type == eImage_uI4)
    return sizeof(int);
  else if (m_eData_type == eImage_realIEEE)
    return sizeof(float);
  else
    return -1;  // For other data types, we don't know the bytes per pixel.
}

int Cimage::nDataSizeNeeded(void)
{
 if ( (0 < m_nDim[0]) && (0 < m_nDim[1]) )
   return (m_nDim[0] * m_nDim[1] * nBytesPerPixel());
 else
   return -1;
}

int Cimage::nGetDataType(void)
{
  return m_eData_type;
}

Cstring Cimage::sGetDataType(void)
{
  if (m_eData_type == eImage_byte)
    return ((Cstring)D_K_SignedChar);
  else if (m_eData_type == eImage_ubyte)
    return ((Cstring)D_K_UnsignedChar);
  else if (m_eData_type == eImage_I2)
    return ((Cstring)D_K_ShortInt);
  else if (m_eData_type == eImage_I4)
    return ((Cstring)D_K_LongInt);
  else if (m_eData_type == eImage_uI2)
    return ((Cstring)D_K_UnsignedShortInt);
  else if (m_eData_type == eImage_uI4)
    return ((Cstring)D_K_UnsignedLongInt);
  else if (m_eData_type == eImage_realIEEE)
    return ((Cstring)D_K_FloatIEEE);
  else if (m_eData_type == eImage_compressed)
    return ((Cstring)D_K_Compressed);
  else
    return ((Cstring)D_K_OtherType);
}

int
Cimage::nAvgSD(const int nPxArea[4],
               const float fExcludeMin,
               const float fExcludeMax,
               float *pfAvg, float *pfSD, float *pfMin, float *pfMax)
{

  // Compute min, max, average and standard deviation in a portion of an image.
  // Exclude from the average calculation pixels with
  // values <= fExcludeMin and >= fExcludeMax (but not from
  // min, max determination.
  // Look only in the the image area defined by nPxArea.

  int   i, j;
  int   nStat;
  float fValue;
  float fSum1, fSum2, fNum;

  fSum1 = 0.0;
  fSum2 = 0.0;
  fNum  = 0.0;

  *pfMin = (this->*prfGetPixel)(nPxArea[0], nPxArea[1]);
  *pfMax = (this->*prfGetPixel)(nPxArea[0], nPxArea[1]);
  for (j = nPxArea[1]; j < nPxArea[1]+nPxArea[3]; j++)
    {
      for (i = nPxArea[0]; i < nPxArea[0]+nPxArea[2]; i++)
        {
          fValue = (this->*prfGetPixel)(i, j);
          if ( (fValue >= fExcludeMin) && (fValue <= fExcludeMax) )
            {
              // This pixel passes all tests to include in the background and sd
              //    calculation.
              fSum1 = fSum1 + fValue;
              fSum2 = fSum2 + (fValue * fValue);
              fNum  = fNum + 1.0f;
            }
          if (fValue < *pfMin) *pfMin = fValue;
          if (fValue > *pfMax) *pfMax = fValue;
        } // end of i loop
    }  // end of j loop

  // For a valid average and sd to be calculated we must have at least 1
  // pixel in the average and more than 1/2 of the pixels in the specified
  // area must contribute!

  if ( (fNum > 1.0) && (fNum > (nPxArea[2] * nPxArea[3] / 2)) )
    {
      *pfAvg = fSum1 / fNum;
      fSum2  = (fSum2 - (fSum1*fSum1 / fNum)) / (fNum - 1.0f);
      if (fSum2 >= 0.0)
        {
          *pfSD = (float)sqrt((double)fSum2);
          nStat = 0;
        }
      else
        {
          *pfSD = (float)sqrt((double)-fSum2);
//        cout << "nAvgSD, sqrt of negative!\n";
          nStat = 1;
        }
    }
  else
    {
      //    cout << "nAvgSD, not enough pixels to compute!\n";
      nStat = 2;
    }
  return (nStat);
}

float
Cimage::fCvtUItoFloat(const unsigned short int uiValue)
{
  return ((float)uiValue);
}

float
Cimage::fCvtRAXIStoFloat(const unsigned short int uiValue)
{
  // Convert an RAXIS style pixel value passed in as unsigned short int
  // to a floating point number.  The pixel value is not necessarily
  // part of this image's data.

  if (32767 >= uiValue)
     return ((float) uiValue);
  else
    {
      return (((float)((short int)(uiValue))
               + m_fCompressInfo[1]) * m_fCompressInfo[0]);
    }
}

float
Cimage::fCvtCBFtoFloat(const unsigned short int uiValue)
{
  // Convert a CBF style pixel value passed in as unsigned short int
  // to a floating point number.  The pixel value is not necessarily
  // part of this image's data.

  return (  (float)((short int)(uiValue)));
}

float
Cimage::fCvtCBFtoFloat(const long int nValue)
{
  // Convert a CBF style pixel value passed in as long int
  // to a floating point number.  The pixel value is not necessarily
  // part of this image's data.

  return (  (float)nValue );
}

float
Cimage::fGetRAXISPixel(const int nWhere0, const int nWhere1)
{
  // Get a pixel at location nWhere[0], nWhere[1]
  // Pixel values above 32767 are compressed by a factor of m_fCompressInfo[0].
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (nWhere0 < 0)
      || (nWhere1 < 0)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (-1.0);
#endif

  // Is there any advantage to this (int) casting?
  // Maybe cast directly to float?

  int nTemp;
  nTemp = (int) *(m_The_Data.pi+nWhere1*m_nDim[0] + nWhere0);
  if (0 <= nTemp)
     return ((float) nTemp);
  else
     return (((float)(nTemp) + m_fCompressInfo[1]) * m_fCompressInfo[0]);
}

float
Cimage::fGetCBFPixel(const int nWhere0, const int nWhere1)
{
  // Get a pixel at location nWhere[0], nWhere[1]
  // Pixel values above 32767 are compressed by a factor of m_fCompressInfo[0].
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (nWhere0 < 0)
      || (nWhere1 < 0)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (-1.0);
#endif

  // Is there any advantage to this (int) casting?
  // Maybe cast directly to float?

  int nTemp;
  nTemp = (int) *(m_The_Data.pn+nWhere1*m_nDim[0] + nWhere0);
  return ((float) (nTemp - m_nPixValOffset));
}


float
Cimage::fCvtDIP2030toFloat(const unsigned short int uiValue)
{
  // Convert a DIP2030 style pixel value passed in as unsigned short int
  // to a floating point number.  The pixel value is not necessarily
  // part of this image's data.

  if (32767 >= uiValue)
     return ((float) uiValue);
  else
    {
      if (5.0f == m_fCompressInfo[0])
        {
          return ( (float)( (int)(~((short int)uiValue)) << 5));
          //      return ( (float)(32 * ( -1 - (short int)(uiValue))));
        }
      else
        {
          return ( (float)(( (int)(~((short int)uiValue)) << 8)+32768));
        }
    }
}

float
Cimage::fGetDIP2030Pixel(const int nWhere0, const int nWhere1)
{
  // Get a pixel at location nWhere[0], nWhere[1]
  // Pixel values above 32767 are compressed by a factor of m_fCompressInfo[0].
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (nWhere0 < 0)
      || (nWhere1 < 0)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (-1.0);
#endif

  // Is there any advantage to this (int) casting?
  // Maybe cast directly to float?

  short int iTemp;
  iTemp =  *(m_The_Data.pi+nWhere1*m_nDim[0] + nWhere0);
  if (0 <= iTemp)
     return ((float) iTemp);
  else
    {
      if (5.0f == m_fCompressInfo[0])
        {
          return ( (float)( (int)(~((short int)iTemp)) << 5));
          //  return ( (float)(32 * (-1 - iTemp)));
        }
      else
        {
          return ( (float)((int)(~(iTemp) << 8)+32768));
        }
    }
}

float
Cimage::fGetSiemensPixel(const int nWhere0, const int nWhere1)
{
  // Get a pixel at location nWhere[0], nWhere[1]
  // Pixel values equal to 65535 may be in overflow table.
#ifndef DISABLE_IMAGE_BOUNDS_CHECKING
  if (   (nWhere0 < 0)
      || (nWhere1 < 0)
      || (nWhere0 >= m_nDim[0])
      || (nWhere1 >= m_nDim[1]) )
    return (-1.0);
#endif

  unsigned short int uiTemp;
  register int nOffset;
  //float fTemp;
  nOffset = nWhere1 * m_nDim[0] + nWhere0;
  uiTemp =  *(m_The_Data.pui + nOffset);

  if (65535 == uiTemp)
    {
      register int i;
      for (i = 0; i < m_nOverflowCount; i++)
        {
          if (m_pnOverflowOffset[i] == nOffset)
            {
              return (m_pfOverflowValue[i]);
            }
        }
      return ((float)uiTemp);
    }
  else
    {
      return ((float)uiTemp);
    }
}

float
Cimage::fCvtSiemenstoFloat(const unsigned short int uiValue)
{
  // BOGUS! NOT WORKING YET!
  // Convert a Siemens style pixel value passed in as unsigned short int
  // to a floating point number.  The pixel value is not necessarily
  // part of this image's data.

  if (65535 > uiValue)
     return ((float) uiValue);
  else
    {
     return ((float) uiValue);
    }
}

int
Cimage::nSetHeader(Cimage_header *poHeader)
{
  // Change the header of the image completely.
  // Update Cimage:: member variables as needed (data type, satval, etc)
  // Reallocate memory if needed, otherwise data remains intact.

  int     nStat;
  Cstring sTemp;
  float   fTemp;

  nStat = 0;

  if (NULL != poHeader)
    {
      // Copy new header to this objects member var header

      m_oHeader = *poHeader;
    }

  // Interpret the header

  (void) m_oHeader.nGetValue(Cimage_header::ms_sSize1, &m_nDim[0]);
  (void) m_oHeader.nGetValue(Cimage_header::ms_sSize2, &m_nDim[1]);

  if (1 < ms_nVerbose)
    {
      cout << "     Image is " << m_nDim[0] << " by " << m_nDim[1] << " pixels.\n";
    }

  if (0 == m_oHeader.nGetValue(Cimage_header::ms_sDataType, &sTemp))
    {
      if (D_K_ShortInt == sTemp)
        {
          m_eData_type = eImage_I2;
          m_fSatVal  = 32767.0;
        }
      else if (D_K_UnsignedShortInt == sTemp)
        {
          m_eData_type  = eImage_uI2;
          prfGetPixel = &Cimage::fGetuiPixel;
          m_fSatVal   = 65535.0;
        }
      else if (D_K_LongInt == sTemp)
        {
          m_eData_type = eImage_I4;
          prfGetPixel        = &Cimage::fGetCBFPixel;
          prfCvtUItoFloat    = &Cimage::fCvtCBFtoFloat;
          prnSetPixel        = &Cimage::nSetCBFPixel;
          m_nPixValOffset    = ms_nPilatusPixelValueOffset;
          m_fSatVal   = (float)(2147483648.0 * 2.0 - 1.0);
        }
      else if (D_K_UnsignedLongInt == sTemp)
        {
          m_eData_type = eImage_uI4;
          m_fSatVal   = (float)(2147483648.0 - 1.0);
        }
      else if (D_K_SignedChar == sTemp)
        {
          m_eData_type = eImage_byte;
          m_fSatVal   = 127.0;
        }
      else if (D_K_UnsignedChar == sTemp)
        {
          m_eData_type = eImage_ubyte;
          m_fSatVal   = 255.0;
        }
      else if ( (D_K_FloatIEEE == sTemp) || ("float" == sTemp) )
        {
          m_eData_type = eImage_realIEEE;
          prfGetPixel = &Cimage::fGetfPixel;    
          m_fSatVal = 4194176.0;  // This could be incorrect
        }
      else if (D_K_Compressed == sTemp)
        {
          m_eData_type = eImage_compressed;
          m_fSatVal   = 65535.0;
        }
      else
        {
          m_eData_type = eImage_other_type;
          m_fSatVal   = 65535.0;
        }
    }
  else if (0 == m_oHeader.nGetValue(Cimage_header::ms_sType, &sTemp))
    {
      // Allow some backwards compatible stuff from Marty Stanton

      if (   (D_K_MadMarty == sTemp)
          || (D_K_UnsignedShortIntMarty == sTemp) )
        {
          // Marty's mad is really unsigned short int

          m_eData_type = eImage_uI2;
          prfGetPixel = &Cimage::fGetuiPixel;
          m_fSatVal   = 65535.0;
        }
       else if ( (D_K_LongIntMarty == sTemp) || (D_K_SignedLongMarty == sTemp) )
        {
          m_eData_type = eImage_I4;
          prfGetPixel        = &Cimage::fGetCBFPixel;
          prfCvtUItoFloat    = &Cimage::fCvtCBFtoFloat;
          prnSetPixel        = &Cimage::nSetCBFPixel;
          m_nPixValOffset    = ms_nPilatusPixelValueOffset;
          m_fSatVal   = (float)(2147483648.0 * 2.0 - 1.0);
        }
      else if ( (D_K_FloatIEEE == sTemp) || ("float" == sTemp) )
        {
	  // Allow some backwards compatible stuff from Marty Stanton
          m_eData_type = eImage_realIEEE;
          prfGetPixel = &Cimage::fGetfPixel;    
          m_fSatVal = 4194176.0;  // This could be incorrect
        }
      else
        {
          m_eData_type = eImage_other_type;
          m_fSatVal   = 65535.0;
        }
    }
  else
    {
      // Default data_type is unsigned short int

      m_eData_type  = eImage_uI2;
      prfGetPixel = &Cimage::fGetuiPixel;
      m_fSatVal   = 65535.0;
//      cout << "     Data_type defaults to " << sGetDataType() << ".\n";
    }

  if (1 < ms_nVerbose)
    {
      cout << "     Data_type is " << sGetDataType() << ".\n";
    }

  if (0 == m_oHeader.nGetValue(Cimage_header::ms_sCompression, &sTemp))
    {
      if (1 < ms_nVerbose)
        {
          cout << "     Compression_type is " << sTemp << ".\n";
        }
      if ("MSC1" == sTemp)
      {
          m_eCompression = eImage_MSC1_compression;
      }
      else if ("None" == sTemp)
      {
          m_eCompression = eImage_no_compression;
      }
      else if ("PCK" == sTemp)
      {
          m_eCompression = eImage_PCK1_compression;
      }
      else if ("WinBMP" == sTemp)
      {
          m_eCompression = eImage_other_compression;
      }
#ifdef RAXIS
      else if ("RAXIS" == sTemp)
        {
          m_eCompression = eImage_MSC1_compression;
          (void) m_oHeader.nGetValue( Cstring ("RAXIS_DETNUM"), &m_nDetNum);
        }
#endif
      else if ("PCK" == sTemp)
      {
          m_eCompression = eImage_PCK1_compression;
      }
#ifdef SIEMENS
      else if ("BS1" == sTemp)
      {
          m_eCompression = eImage_BS1_compression;
      }
#endif
      else if ("CBF" == sTemp)
      {
          m_eCompression = eImage_CBF_compression;
      }
#ifdef BAS2000
      else if ("BAS2000" == sTemp)
      {
          m_eCompression = eImage_BAS2000_compression;
      }
#endif
#ifdef HIPIC
      else if ("HIPIC" == sTemp)
      {
          m_eCompression = eImage_HIPIC_compression;
      }
#endif
      else
      {
          m_eCompression = eImage_other_compression;
      }
    }
  else
    {
      m_eCompression = eImage_no_compression;
    }

  if (0 == m_oHeader.nGetValue(Cimage_header::ms_sByteOrder, &sTemp))
    {
      if (1 < ms_nVerbose)
        {
          cout << "     Byte_order is " << sTemp << ".\n";
        }

      if (D_K_BigEndian == sTemp)
        m_eImgByteOrder = eCPU_big_endian;
      else if (D_K_LittleEndian == sTemp)
        m_eImgByteOrder = eCPU_little_endian;
      else
        m_eImgByteOrder = eCPU_other_endian;
    }
  else
    {
      m_eImgByteOrder = (eCPU_byte_orders)nGetByteOrderCPU();
      if (1 < ms_nVerbose)
        {
          cout << "     No BYTE_ORDER keyword in image header!\n"
               << "        Byte order is assumed to be the byte order of this CPU.\n";
        }
    }

  // Even though the default saturated value is set above based on data type,
  // check if the header has a definition for the saturated value here

  if (0 == m_oHeader.nGetValue(Cimage_header::ms_sSaturatedValue, &fTemp))
    {
      m_fSatVal = fTemp;
    }

#ifdef RAXIS
  // Allow a standard d*TREK image to use R-AXIS style pixel
  // compression.  This should be moved below to apply to all
  // image types.
//+2013-05-12 JWP
// DANGER, if RAXIS_COMPRESSION_RATIO is not 8, 32, 128 (e.g. if it is 1)
// or if the data_type is not unsigned short int, then this does not apply

  if (eImage_uI2 == nGetDataType())
    {
      if (0 == m_oHeader.nGetValue(Cimage_header::ms_sRaxisCompressionRatio,
			       &m_fCompressInfo[0])) 
	{
//-2013-05-12 JWP
      if ((float)1.0 >= m_fCompressInfo[0])
        {
          m_fCompressInfo[1] = 0.0;
        }
      else
        {
          m_fCompressInfo[1] = 32768.0;
        }
      nStat = m_oHeader.nGetValue(D_K_SaturatedValue, &m_fSatVal);
      if (0 != nStat)
        {
          // Saturated value is not in the image header so force a value
          nStat = 0;
          if ((float)8.0 == m_fCompressInfo[0])
            // It is an raxis2
            m_fSatVal = 262136.0;                 // 2**18 - 2**3
          else if ((float)32.0 == m_fCompressInfo[0])
            // It is an raxis4
            //m_fSatVal = 1048544.0;               //  2**20 - 2**5
            // Per Katsu Sasaki
            m_fSatVal = 1048512.0;               //  2**20 - 2**5 - 2**5
          else if ((float)128.0 == m_fCompressInfo[0])
            m_fSatVal = 4194176.0;               //  2**22 - 2**7
        }
      prfGetPixel      = &Cimage::fGetRAXISPixel;
      prfCvtUItoFloat  = &Cimage::fCvtRAXIStoFloat;

      // RB
      //if( eImage_no_compression == m_eCompression )
      //  m_eCompression = eImage_MSC1_compression;
      /// end RB
    }
    }

#endif
/*
   set m_sfilename
*/

  if (NULL != poHeader)
    {
      // If a new header was specified, check that enough space is
      // allocated for the pixels.  If not, go get new space

      int nBytes = nDataSizeNeeded();
      if (m_nBytesAllocated < nBytes)
        {
          if (NULL != m_The_Data.pc)
            delete [] m_The_Data.pc;
          m_The_Data.pc = NULL;
          m_The_Data.pc = new char [nBytes]; // pc & pn are the same via union
          m_nBytesAllocated = nBytes;
        }
    }

  return (nStat);
}

// Bin pixel intensities acording to a custom bin array vecBins, and a "test-bit" bin array vecBitBins. 
// The "test-bit" bin array is useful to check if a certain bit is lacking in integer pixel values.
// Bin arrays must start with a bin to count values less than the minimum value of the bin array. 
// Bin arrays must end with a bin to count values equal or greater than the maximum value of the bin array..
void Cimage::vBuildIntensityHistogram(std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins, 
                                      std::vector<int>& vecBitBins, 
                                      int& nNumberOfPixValuesInOverlapRange_PMT1,
                                      int& nNumberOfPixValuesInOverlapRange_PMT2,
                                      DTREK_WORD wCtrl)
{
    int         nPixelCount = m_nDim[0] * m_nDim[1];    
    int         nPixIntensity = 0;

    nNumberOfPixValuesInOverlapRange_PMT1 = 0;
    nNumberOfPixValuesInOverlapRange_PMT2 = 0;
    
    ////////////////////////////////////////////////////////////////////
    int         nBitShiftConstant = -1;
    if( !(wCtrl & DTREK_IMAGE_BIH_CCD) )
    {
        nBitShiftConstant = wCtrl & DTREK_IMAGE_BIH_RAXIS_II ? 8 : 32;
    }
    ////////////////////////////////////////////////////////////////////

    for(int ii=0; ii < m_nDim[0]; ii++)
    {
        for(int jj=0; jj < m_nDim[1]; jj++)
        {
            nPixIntensity = (int)(this->*prfGetPixel)(ii, jj);

            ////////////////////////////////////////////////////////////////////////////////////////////////
            if( !(wCtrl & DTREK_IMAGE_BIH_CCD) )
            {
                if( nPixIntensity >= PMT_OVERLAP_BEGIN && nPixIntensity <= PMT_OVERLAP_END )
                {
                    if( nPixIntensity / nBitShiftConstant * nBitShiftConstant == nPixIntensity ) // diviseable
                        nNumberOfPixValuesInOverlapRange_PMT2++;
                    else
                        nNumberOfPixValuesInOverlapRange_PMT1++;
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            if( wCtrl & DTREK_IMAGE_BIH_CUSTOM && vecBins.size() > 0 )
            {
                if( nPixIntensity < vecBins[0].m_nEnd )
                    vecBins[0].m_nCount += 1;
                else if( nPixIntensity >= vecBins[vecBins.size()-1].m_nBegin )
                {
                    vecBins[vecBins.size()-1].m_nCount += 1;
                }
                else
                {
                    for(int kk=1; kk < vecBins.size()-1; kk++)
                    {
                        if( nPixIntensity >= vecBins[kk].m_nBegin && nPixIntensity < vecBins[kk].m_nEnd )
                        {
                            vecBins[kk].m_nCount += 1;
                            break;
                        }
                    }
                }
            }

            
            if( wCtrl & DTREK_IMAGE_BIH_BITTEST && vecBitBins.size() > 0 )
            {
                bool    bAtLeastOneBitFound = false;

                if( nPixIntensity == 0 )
                    vecBitBins[0] += 1;
                else
                {
                    int     nMask = 0;
                    for(int nBitNumber = 1; nBitNumber < vecBitBins.size()-1; nBitNumber++)
                    {
                        nMask = 1 << (nBitNumber - 1);
                
                        if( nMask & nPixIntensity )
                        {
                            vecBitBins[nBitNumber] += 1;
                            bAtLeastOneBitFound = true;
                        }
                    }
                }
                
                if( nPixIntensity != 0 && !bAtLeastOneBitFound )
                    vecBitBins[vecBitBins.size()-1] += 1;
            }
        }  // jj
    } // ii
}

int 
Cimage::nConvertIntToUIShort(const int nZeroConversion)
{
  // Convert "in place" an image pixel data of signed long to unsigned short int.
  // Typically this will be used by the Cnonunf class to convert a long int image into unsigned
  // short int.  Pixel values less than 0 will be set to 0.  Pixel values of 0 will be set to 
  // nZeroConversion which has a default value of 1.

  // WARNING: Do not use this routine except for converting a Cnonunf mask file to a uishort
  //          without first checking for suitability.

  int nStat;
  unsigned short int *puiNewData;
  unsigned short int uiTemp;
  float fPix;
  unsigned int uiConvertedZero;
  uiConvertedZero = (unsigned short int)nZeroConversion;
  if (1 == m_nPixValOffset)
    uiConvertedZero = 1;
      
  cout << "INFO: 0's converted to a value of " << uiConvertedZero << endl;
  if (eImage_I4 != m_eData_type)
    {
      cout << "WARNING, attempt to convert non-I4 image to uiShort image\n";
      return (-1);
    }

  nStat = 0;
  puiNewData = m_The_Data.pui;
  int i, j;
  for (j = 0; j < m_nDim[1]; j++)
    {
      for (i = 0; i < m_nDim[0]; i++)
        {
          fPix =  (this->*prfGetPixel)(i, j);  // This will use m_nPixValOffset
          if (0.0 > fPix) 
            uiTemp = 0;
          else if (0.0 == fPix) 
            uiTemp = uiConvertedZero;
          else if (65534 < fPix)
            // This image cannot be used for anything except a mask
            uiTemp = 65535;
          else
            uiTemp = (unsigned short int) fPix + m_nPixValOffset;
          *puiNewData++ = uiTemp;
        }
    }

  // Change member variables

  m_eData_type = eImage_uI2;
  m_nPixValOffset = 0;
  m_fSatVal = 65535;

  // Change the member pointer functions

  prfGetPixel      = &Cimage::fGetPixel;
  prfCvtUItoFloat  = &Cimage::fCvtUItoFloat;
  prnSetPixel      = &Cimage::nSetPixel;

  return (nStat);
}

int 
Cimage::nUnbin(void)
{
  int nStat = 0;

  // Probably should do this only if (eImage_uI2 != m_eData_type)

  // Unbin the image data and adjust the image header if needed.
  //cout << "DTREK_UNBIN is set B\n";

  // Check if enough memory is allocated in order to unbin in place
  int nBytes = 0;

  nBytes = nDataSizeNeeded();

  if (m_nBytesAllocated < (4 * nBytes))
    {
      // ERROR
      
      cout << "ERROR in Cimage::nUnbin, not enough memory allocated.\n";
      nStat = 1;
      return (nStat);
    }

  nStat = m_oHeader.nUnbinHeader();

  if (-1 == nStat)
    {
      // Header was already unbinned, so be careful.
    }

  // Now copy the data in place

  float fPixS;
  int   i, j;
  int   nDim0, nDim1;
  nDim0 = m_nDim[0];
  nDim1 = m_nDim[1];
  for (j = nDim1-1; j >= 0; j--)
    {
      for (i = nDim0-1; i >=0; i--)
        {
          // This should work since there are no array bounds checking in 
          // prfGetPixel and prnSetPixel
          
          // Restore m_nDim[0] and m_nDim[1],
          m_nDim[0] = nDim0;
          m_nDim[1] = nDim1;

          fPixS = (this->*this->prfGetPixel)(i, j);
                  
          // m_nDim[0] and m_nDim[1] are used in prnSetPixel

          m_nDim[0] = 2 * m_nDim[0];
          m_nDim[1] = 2 * m_nDim[1];
          (this->*this->prnSetPixel)(2*i,   2*j,   fPixS);
          (this->*this->prnSetPixel)(2*i+1, 2*j,   fPixS);
          (this->*this->prnSetPixel)(2*i,   2*j+1, fPixS);
          (this->*this->prnSetPixel)(2*i+1, 2*j+1, fPixS);
        }
    }
  m_nData_size = m_nBytesAllocated;

  // Finally change the keywords for the array size in the header

  m_nDim[0] = 2 * nDim0;
  m_nDim[1] = 2 * nDim1;

	m_oHeader.nReplaceValue(Cstring(D_K_Size1), m_nDim[0]);
	m_oHeader.nReplaceValue(Cstring(D_K_Size2), m_nDim[1]);

	return (nStat);
}

