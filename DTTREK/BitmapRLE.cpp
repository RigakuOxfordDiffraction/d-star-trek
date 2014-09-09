// 
// Copyright © 2001 Rigaku/MSC, Inc.
//                  9009 New Trails Drive  
//                  The Woodlands, TX, USA  77381  
// 
// The contents are unpublished proprietary source  
// code of RigakuMSC, Inc.  
// 
// All rights reserved   
// 
// BitmapRLE.cpp    Initial author: T.L.Hendrixson           Jul 2001
//
// Description:
//
//    This file contains the implementation of the BitmapRLE class, which
//    provides generation and interpretation of bitmaps stored using a
//    modified run-length encoding (RLE) scheme.
//
// ToDo:
//
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

/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#ifdef WIN32
#pragma warning(disable:4786) // so warning about debug term > 255 chars ignored
#endif

#if (!defined(ICC)) && (!defined(SUNOS))
#include <limits>
#else
#include <limits.h>
#endif

#if (defined(_MSC_VER) && (_MSC_VER==1100))
#include <fstream.h>
#else
#include <fstream>
#endif

#include "BitmapRLE.h"

// Can't include this here because a macro definition somewhere of max() will 
// mess up the std::numberic_limits<T>::max() routine.
//#include "Cimage.h"

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                         Static Member Variables                          *
 ****************************************************************************/

#if (!defined(ICC)) && (!defined(SUNOS))
const unsigned short BitmapRLE::m_MaximumRunLength = (std::numeric_limits<unsigned short>::max()>>1);
const unsigned short BitmapRLE::m_MSB = ((std::numeric_limits<unsigned short>::max()>>1)+1);
#else
const unsigned short BitmapRLE::m_MaximumRunLength = 32767;
const unsigned short BitmapRLE::m_MSB = 32768;
#endif
/****************************************************************************
 *                             Global variables                             *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                           Function prototypes                            *
 ****************************************************************************/

/****************************************************************************
 *                           Non-class functions                            *
 ****************************************************************************/

/****************************************************************************
 *                        Additional Include Files                          *
 ****************************************************************************/

#include "Cimage.h"

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                         Public Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************/
BitmapRLE::BitmapRLE(void) :
   m_BigEndian(true),
   m_pBuffer(NULL),
   m_AllocatedSize(0),
   m_BufferSize(0),
   m_ResultantImageSize(0)
{
   unsigned short one = 1;
   char *p = reinterpret_cast<char *>(&one);
   m_BigEndian = (0 == static_cast<int>(*p));
}

/****************************************************************************
 ****************************************************************************/
BitmapRLE::~BitmapRLE()
{
   delete [] m_pBuffer;
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::Apply(Cimage &Image)
{
   if(eImage_uI2 != Image.m_eData_type && eImage_realIEEE != 
   Image.m_eData_type)
      return false;

   if(Image.nGetDimension() != m_ResultantImageSize)
      return false;

   // Apply the bitmap image to the image passed in.  Only pixels whose
   // bitmap pixel is 0 will be modified.

   Image.nSetNextPixel(0,0);

   // Don't forget that we have to start from after the marker.

   unsigned short *pEnd = m_pBuffer+m_BufferSize;
   for( unsigned short *p = m_pBuffer+2; p < pEnd; p++){

      long n = GetRunLength(*p);

      if(eImage_uI2 == Image.m_eData_type){
         if(!GetValue(*p)){
            // Have to do this to force the correct overloaded
            // function to be called.
            unsigned short zero = 0;
            for(int j = 0; j < n; j++)
               Image.vSetNextPixel(zero);
         }
         else{
            // Easiest way to increment internal image pointer is to
            // get pixel values.  Maybe not the fastest, but the easiest.
            for(int j = 0; j < n; j++)
               Image.uiGetNextPixel();
         }
      }
      else if(eImage_realIEEE == Image.m_eData_type){
         if(!GetValue(*p)){
            // Have to do this to force the correct overloaded
            // function to be called.
            float zero = 0.0f;
            for(int j = 0; j < n; j++)
               Image.vSetNextPixel(zero);
         }
         else{
            // Easiest way to increment internal image pointer is to
            // get pixel values.  Maybe not the fastest, but the easiest.
            for(int j = 0; j < n; j++)
               Image.fGetNextPixel();
         }
      }
   }

   return true;
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::Generate(Cimage &Image)
{

   // NOTE this routine only works for images of type unsigned short
   // and float.

   if(eImage_uI2 != Image.m_eData_type && eImage_realIEEE != 
   Image.m_eData_type)
      return false;

   // In a worst case scenario, we have a alternating series of pixels (
   // (0, x, 0, y, ...) which means the bitmap is the same number of
   // characters as there are pixels in the image.

   long size = Image.nGetDimension();

   unsigned short *pBuffer = new unsigned short [(size+2)*sizeof(unsigned short)];

   // Fill the bitmap.

   unsigned short *pBufferEnd;
   
   // We can't use a template for different image types because of the
   // design of the Cimage class, so use separate routines.

   if(eImage_uI2 == Image.m_eData_type)
      pBufferEnd = GenerateFromUnsignedShortImage(Image,pBuffer);
   else if(eImage_realIEEE == Image.m_eData_type)
      pBufferEnd = GenerateFromFloatImage(Image,pBuffer);

   size = pBufferEnd-pBuffer;
   Set(pBuffer,size);

   delete pBuffer;

   return true;
}

/****************************************************************************
 ****************************************************************************/
void BitmapRLE::Get(unsigned short * const pBuffer)
{
   if(NULL != pBuffer){
      memcpy(reinterpret_cast<void *>(pBuffer),
             reinterpret_cast<void *>(m_pBuffer),
             m_BufferSize*sizeof(unsigned short));
   }
   return;
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::Interpret(Cimage &Image)
{
   if(eImage_uI2 != Image.m_eData_type && eImage_realIEEE != 
   Image.m_eData_type)
      return false;

   if(Image.nGetDimension() != m_ResultantImageSize)
      return false;

   // Create the bitmap image.

   Image.nSetNextPixel(0,0);

   // Don't forget that we have to start from after the marker.

   unsigned short *pEnd = m_pBuffer+m_BufferSize;
   for(unsigned short *p = m_pBuffer+2; p < pEnd; p++){

      long n = GetRunLength(*p);

      if(eImage_uI2 == Image.m_eData_type){
         unsigned short value = (GetValue(*p)) ? 1 : 0;
         for(int j = 0; j < n; j++)
            Image.vSetNextPixel(value);
      }
      else if(eImage_realIEEE == Image.m_eData_type){
         float value = (GetValue(*p)) ? 1.0 : 0.0;
         for(int j = 0; j < n; j++)
            Image.vSetNextPixel(value);
      }
   }

   return true;
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::Set(unsigned short * const pBuffer,
                    const long Size)
{
   if(NULL == pBuffer || 0 >= Size)
      return false;

   // Make sure that the buffer input has the correct marker.

   char *p = reinterpret_cast<char *>(pBuffer);
   if('B' != p[0] || 'R' != p[1] || 'L' != p[2] || 'E' != p[3])
      return false;

   // Make sure we have enough space to hold the bitmap.

   if(m_AllocatedSize < Size){
      delete [] m_pBuffer;
      m_AllocatedSize = Size;
      m_pBuffer = new unsigned short [m_AllocatedSize];
   }
   m_BufferSize = Size;

   memcpy(reinterpret_cast<void *>(m_pBuffer),
          reinterpret_cast<void *>(pBuffer),
          m_BufferSize*sizeof(unsigned short));

   m_ResultantImageSize = 0;
   unsigned short *pEnd = pBuffer+Size;
   for(unsigned short *ptr = pBuffer+2; ptr < pEnd; ptr++)
      m_ResultantImageSize += GetRunLength(*ptr);

   return true;
}

/****************************************************************************
 ****************************************************************************/
long BitmapRLE::Size(void)
{
   // Don't forget to account for the marker.
   return m_BufferSize-2;
}

/****************************************************************************
 ****************************************************************************/
long BitmapRLE::SizeOf(void)
{
   return m_BufferSize;
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::Write(const std::string &Filename)
{

#if (defined(_MSC_VER) && (_MSC_VER==1100))
   ofstream out(Filename.c_str(),ios::out|ios::app|
                     ios::binary);

   if(out.good())
      out.write(const_cast<const char *>(reinterpret_cast<char *>(m_pBuffer)),
                m_BufferSize*sizeof(unsigned short));

#else

#if (!defined(ICC)) && (!defined(SUNOS))
   std::ofstream out(Filename.c_str(),std::ios_base::out|std::ios_base::app|
                     std::ios_base::binary);
#else
   std::ofstream out(Filename.c_str(),std::ios::out|std::ios::app|
                     std::ios::binary);
#endif

   if(out.good())
      out.write(const_cast<const char *>(reinterpret_cast<char *>(m_pBuffer)),
                m_BufferSize*sizeof(unsigned short));
#endif

   bool OK = out.good();

   if(out.is_open())
      out.close();

   return OK;
}

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                       Protected Member Functions                       **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/
  
/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                        Private Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************/
void BitmapRLE::AddMarker(unsigned short * const pBuffer)
{
   // Put the marker in the compression buffer.

   char *pMarker = reinterpret_cast<char *>(pBuffer);
   pMarker[0] = 'B';
   pMarker[1] = 'R';
   pMarker[2] = 'L';
   pMarker[3] = 'E';

   return;
}

/****************************************************************************
 ****************************************************************************/
unsigned short *BitmapRLE::GenerateFromFloatImage(Cimage &Image,
                                                  unsigned short * const pBuffer)
{
   // Start filling up the bitmap.

   AddMarker(pBuffer);

   float *pData = reinterpret_cast<float *>(Image.pvGetData());
   float *pEnd = pData+Image.nGetDimension();
  
   // Init the variables

   unsigned short *p = pBuffer+2;
   long cnt = 1;
   bool ones = (0.0 != *pData);
   pData++;

   for(NULL; pData < pEnd; pData++){

      // If current pixel is "same" as the previous pixel, just increment
      // the count.
      if((ones && (0.0 != *pData)) || (!ones && (0.0 == *pData)))
         cnt++;
      // Otherwise, the pixel value has switched, so we will add our 
      // current info to the internal buffer
      else{
         p = SetRunLengthValue(p,cnt,ones);
         cnt = 1;
         ones = (0.0 != *pData);
      }
   }

   // Now that we are done, don't forget to put our last count in the
   // buffer.

   return SetRunLengthValue(p,cnt,ones);
}

/****************************************************************************
 ****************************************************************************/
unsigned short *BitmapRLE::GenerateFromUnsignedShortImage(Cimage &Image,
                                                          unsigned short * const pBuffer)
{
   // Start filling up the bitmap.

   AddMarker(pBuffer);

   unsigned short *pData = reinterpret_cast<unsigned short *>(Image.pvGetData());
   unsigned short *pEnd = pData+Image.nGetDimension();
  
   // Init the variables

   unsigned short *p = pBuffer+2;
   long cnt = 1;
   bool ones = (0 != *pData);
   pData++;

   for(NULL; pData < pEnd; pData++){

      // If current pixel is "same" as the previous pixel, just increment
      // the count.
      if((ones && (0 != *pData)) || (!ones && (0 == *pData)))
         cnt++;
      // Otherwise, the pixel value has switched, so we will add our 
      // current info to the internal buffer
      else{
         p = SetRunLengthValue(p,cnt,ones);
         cnt = 1;
         ones = (0 != *pData);
      }
   }

   // Now that we are done, don't forget to put our last count in the
   // buffer.

   return SetRunLengthValue(p,cnt,ones);
}

/****************************************************************************
 ****************************************************************************/
long BitmapRLE::GetRunLength(const unsigned short Data)
{
   unsigned short data = (m_BigEndian) ? Data : Swap(Data);
   return static_cast<long>(m_MaximumRunLength&data);
}

/****************************************************************************
 ****************************************************************************/
bool BitmapRLE::GetValue(const unsigned short Data)
{
   unsigned short data = (m_BigEndian) ? Data : Swap(Data);
   return (0 != (m_MSB&data));
}

/****************************************************************************
 ****************************************************************************/
unsigned short *BitmapRLE::SetRunLengthValue(const unsigned short * const pData,
                                             const long RunLength,
                                             const bool NonZeroValue)
{
   unsigned short *p = const_cast<unsigned short *>(pData);
   long cnt = RunLength;

   for(NULL; cnt > 0; cnt -= m_MaximumRunLength){
      unsigned short value = (cnt > m_MaximumRunLength) ? m_MaximumRunLength :
                                                          cnt;
      if(NonZeroValue)
         value = static_cast<unsigned short>(value)|m_MSB;
      else
         value = static_cast<unsigned short>(value);
      *p++ = (m_BigEndian) ? value : Swap(value);
   }

   return p;
}

/****************************************************************************
 ****************************************************************************/
unsigned short BitmapRLE::Swap(const unsigned short Value)
{
   unsigned short swap = Value;
   char *p = reinterpret_cast<char *>(&swap);

   char tmp = p[0];
   p[0] = p[1];
   p[1] = tmp;

   return swap;
}

