#ifndef BITMAPRLE_H
#define BITMAPRLE_H

//# 
//# Copyright © 2001 Rigaku/MSC, Inc.
//#                  9009 New Trails Drive  
//#                  The Woodlands, TX, USA  77381  
//# 
//# The contents are unpublished proprietary source  
//# code of RigakuMSC, Inc.  
//# 
//# All rights reserved   
//# 
//# BitmapRLE.h      Initial author: T.L.Hendrixson           Jul 2001
//#
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

#include <string>

#include "dtrekdefs.h"

/****************************************************************************
 *                           Forward declarations                           *
 ****************************************************************************/

class Cimage;

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                             BitmapRLE class                              *
 ****************************************************************************/

// <P>
//    The BitmapRLE class is used for generating and interpreting a 
//    bitmap constructed from d<SUP>*</SUP>TREK image.  The bitmap is
//    encoded using a variant on the run-length encoding (RLE) scheme,
//    and uses the simplifying assumption that the original data is
//    either zero or not zero.  Thus, this is a lossy compression
//    scheme, in terms of the original data.
// </P>
// <P>
//    The traditional RLE scheme compresses data by using a repeated 2
//    byte sequence.  One byte contains the character to be operated on
//    and the other byte contains the number of instances, or length
//    of the run.  Thus the character sequence <KBD>AAABCCCCCCDDDDDDDD</KBD>
//    would compress to <KBD>'A' 0x03 'B' 0x01 'C' 0x06 'D' 0x08</KBD>,
//    resulting in an 8 character sequence instead of the original 18 
//    character sequence.
// </P>
// <P>
//    The modified RLE scheme used by this class takes advantage of the
//    fact that the only information that needs to be retained is whether
//    a data point (<I>i.e.</I>, a pixel value) is either zero or not zero.
//    An encoded value contains 2 types of information: 1) whether the value
//    is zero or not zero, and 2) the length of the run.  The most
//    significant bit is used to indicate if the data is zero or not zero.
//    If the bit is set, then the data point is non-zero.  If the bit is
//    not set, then the data point is zero.  The remainder of the encoded
//    value contains the length of the run.  Note that the maximum length of 
//    a run encoded in an encode value is 32767 characters.  For runs of 
//    longer than 32767 data points, multiple encoded values must be used.  
//    Thus a character sequence of <KBD>0x00 0x00 0x00 0x34 0x45 0x43 0x02 
//    0x03</KBD> will be compressed as <KBD>0x0003 0x8005</KBD>.  Note that
//    the encoded values are stored in a big-endian fashion.
// </P>
// <P>
//    The first 4 bytes in the encoded bitmap are the characters 
//    <KBD>BRLE</KBD>, which is a marker indicating that the encoded bitmap
//    uses the Bitmap RLE encoding scheme.  These 4 characters are
//    considered to be part of the encoded bitmap and must be present in 
//    any encoded bitmap given to the class, via the <CODE>Set()</CODE>
//    routine, and will be contained in the encoded bitmap returned by
//    the <CODE>Get()</CODE> routine.
// </P>
// <P>
//    <STRONG>NOTE</STRONG>: In the current implementation of this class,
//    only d<SUP>*</SUP>TREK images of type <KBD>eImage_uI2</KBD> and
//    <KBD>eImage_realIEEE</KBD> (unsigned short and float) are 
//    supported.
// </P>

class DTREK_EXPORT BitmapRLE
{

public:

//# Constructor

      // Create an instance of the class.
   BitmapRLE(void);

//# Destructor

      // The destructor for this class.  This will free any memory that
      // was allocated by the class.
   virtual ~BitmapRLE();

      // This routine applys the run-length encoded bitmap to the image
      // <VAR>Image</VAR>.  When the bitmap is applied, only those pixels 
      // whose corresponding bitmap image value are 0 will be modified,
      // having their values set to zero.  Pixels whose corresponding 
      // bitmap image value are non-zero will not be modified.  If the 
      // value returned is <I>true</I>, then the bitmap was interpreted.
      // If the value returned is <I>false</I>, then the bitmap could not 
      // be interpreted or the image was not of the correct type or of 
      // sufficient size to hold the resultant bitmap.
   bool Apply(Cimage &Image);

      // This routine generates a run-length encoded bitmap of the 
      // d<SUP>*</SUP>TREK image <VAR>Image</VAR> and returns a boolean
      // indicating if the bitmap could be generated.  If the value
      // returned is <I>true</I>, then the bitmap was generated.  If the 
      // value returned is <I>false</I>, then the bitmap could not be
      // generated.
   bool Generate(Cimage &Image);

      // This routine places the run-length encoded bitmap in the memory
      // location pointed to by <VAR>pBuffer</VAR>.  If <VAR>pBuffer</VAR>
      // is <KBD>NULL</KBD>, then this routine is a non-operation.  It is
      // assumed that <VAR>pBuffer</VAR> is large enough to hold the bitmap.
      // If it is not, the results of this routine are undefined.
   void Get(unsigned short * const pBuffer);

      // This routine interprets the run-length encoded bitmap, places
      // the resulting image in the d<SUP>*</SUP>TREK image <VAR>Image</VAR>,
      // and returns a boolean indicating if the bitmap could be interpreted.  
      // If the value returned is <I>true</I>, then the bitmap was 
      // interpreted.  If the value returned is <I>false</I>, then the 
      // bitmap could not be interpreted or the image was not of the
      // correct type or of sufficient size to hold the resultant bitmap.
   bool Interpret(Cimage &Image);

      // <P>
      //   This routine sets the run-length encoded bitmap to be the bitmap
      //   stored in the memory location pointed to by <VAR>pBuffer</VAR>.
      //   <VAR>Size</VAR> is the size of the bitmap, in unsigned shorts.
      //   If the value returned is <I>true</I>, then the bitmap was set
      //   correctly.  If the value returned is <I>false</I>, then the 
      //   bitmap was not set correctly.
      // </P>
      // <P>
      //   If the bitmap is not set correctly, the likely causes are a
      //   NULL value for <VAR>pBuffer</VAR>, a non-positive value for
      //   <VAR>Size</VAR>, or an incorrect (or no) marker at the start
      //   of the bitmap.
      // </P>
      // <P>
      //    If the bitmap is set correctly, the class creates an internal
      //    copy of the bitmap, and so the calling routine is free to
      //    deallocate the memory used for storing the bitmap.
      // </P>
   bool Set(unsigned short * const pBuffer, const long Size);

      // This routine returns the number of encoded values in the current 
      // run-length encoded bitmap.  Note that this value should not be
      // used for computing the size of the memory block required for
      // storing an encoded bitmap, as it does not account for the 
      // marker; use <CODE>SizeOf()</CODE> instead.
   long Size(void);

      // This routine returns the size, in unsigned shorts, of the current
      // run-length encoded bitmap (including the marker).
   long SizeOf(void);

      // This routine writes the run-length encoded bitmap to the file
      // specified by <VAR>Filename</VAR> and returns a boolean indicating
      // if the bitmap could be written successfully.  If the value returned
      // is <I>true</I>, then the bitmap was written successfully.  If the
      // value returned is <I>false</I>, then an error occurred while 
      // writing the bitmap to the file.  If the file already exists,
      // the bitmap will be appended to the end of the file.
	   // <HR>
   bool Write(const std::string &Filename);

protected:

private:

//# There is no copy constructor
      // A non-existant copy constructor.  This is defined as a private
      // member function so that user code cannot contain the following:
      // <PRE><CODE>
      // BitmapRLE instance1;
      // BitmapRLE instance2(instance1);
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
   BitmapRLE(const BitmapRLE &Copy);

//# There is no assignment operator
      // A non-existant assignment operator.  This is defined as a private
      // member function so that user code cannot contain the following:
      // <PRE><CODE>
      // BitmapRLE instance1,instance2;
      // instance2 = instance1;
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
   BitmapRLE &operator=(const BitmapRLE &AssignEqual);

      // This routine adds the marker to the memory location indicated
      // by <VAR>pBuffer</VAR>.  If <VAR>pBuffer</VAR> is NULL, then the
      // results of this routine are undefined.
   void AddMarker(unsigned short * const pBuffer);

      // This routine genenerates the bitmap from the d<SUP>*</SUP>TREK 
      // image <VAR>Image</VAR> of type unsigned short (or 
      // <CODE>eImage_uI2</CODE> in d<SUP>*</SUP>TREK parlance) and 
      // places it starting at the memory location pointed to by 
      // <VAR>pBuffer</VAR>.
   unsigned short *GenerateFromUnsignedShortImage(Cimage &Image,
                                                  unsigned short * const pBuffer);

      // This routine genenerates the bitmap from the d<SUP>*</SUP>TREK 
      // image <VAR>Image</VAR> of type float (or 
      // <CODE>eImage_realIEEE</CODE> in d<SUP>*</SUP>TREK parlance) and 
      // places it starting at the memory location pointed to by 
      // <VAR>pBuffer</VAR>.
   unsigned short *GenerateFromFloatImage(Cimage &Image,
                                          unsigned short * const pBuffer);

      // This routine returns the length of the run for the encoded
      // value <VAR>Data</VAR>.
   long GetRunLength(const unsigned short Data);

      // This routine returns the value of the elements in the run
      // for the encoded value <VAR>Data</VAR>.  If the value returned
      // is <I>true<I>, then the elements are non-zero.  If the value
      // returned is <I>false</I>, then the elements are zero.
   bool GetValue(const unsigned short Data);

      // This routine sets the run length encoded values, starting at
      // <VAR>pData</VAR>, for a run of <VAR>RunLength</VAR> values.  If
      // the value of <VAR>NonZeroValue</VAR> is <I>true</I>, then the
      // run contains non-zero values.  If the value of
      // <VAR>NonZeroValue</VAR> is <I>false</I>, then the run contains
      // zero values.  The value returned is a pointer to the next
      // (unused) position in memory after the encoded values have been
      // set.  If <VAR>pData</VAR> is NULL, then the results of this
      // routine are undefined.
   unsigned short *SetRunLengthValue(const unsigned short * const pData,
                                     const long RunLength,
                                     const bool NonZeroValue);

      // This routine swaps the bytes in <VAR>Value</VAR>, and returns
      // the swapped value.  This converts from big-endian to little-endian,
      // and <I>vice versa</I>.
   unsigned short Swap(const unsigned short Value);

//# Member variables

   bool m_BigEndian;

   unsigned short *m_pBuffer;

   const static unsigned short m_MaximumRunLength;
   const static unsigned short m_MSB;

   long m_AllocatedSize;
   long m_BufferSize;
   long m_ResultantImageSize;

};

#endif /* BITMAPRLE_H */
