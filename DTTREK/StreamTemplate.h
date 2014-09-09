#ifndef STREAMTEMPLATE_H
#define STREAMTEMPLATE_H



//# 
//# Copyright © 2003 Rigaku/MSC, Inc.
//#                  9009 New Trails Drive
//#                  The Woodlands, TX, USA  77381
//# 
//# The contents are unpublished proprietary source
//# code of Rigaku/MSC, Inc.
//# 
//# All rights reserved
//# 
//# StreamTemplate.h     Initial author: T.L.Hendrixson             Jan 2003
//#
//# Description:
//#
//# This file contains templates for using Cstrings with the standard
//# iostreams.
//#
//# ToDo:
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
 *                               Include Files                              *
 ****************************************************************************/

#include <iosfwd>

/****************************************************************************
 *                           Forward Declarations                           *
 ****************************************************************************/

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
 *                           Function prototypes                            *
 ****************************************************************************/
/****************************************************************************
 * This routine inserts a Cstring into a basic_ostream.  Since stringstream
 * is derived from basic_ostream and uses the basic_ostream's << operator,
 * we have to use this to write to stringstreams.
 ****************************************************************************/
template <class A, class B>
inline std::basic_ostream<A,B> &operator<<(std::basic_ostream<A,B> &a, 
                                           const Cstring &b)
{
   a << b.string();
   return a;
}


#endif /* STREAMTEMPLATE_H */
