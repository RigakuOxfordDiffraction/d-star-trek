#ifndef DT_IMAGE_POOL_H
#define DT_IMAGE_POOL_H
//
// Copyright (c) 2004 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// ImagePool.h    Initial author: R.B.   23-June-2004
// This file is the header file for class CImagePool
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
#include "Dtrek.h"
#include "dtreksys.h"
#include "Cimage.h"

//+Definitions and constants
//+Forward class declarations

//+Code begin
class DTREK_EXPORT CImagePool
{
public:
    CImagePool();
    ~CImagePool();

    static CImagePool*  poGetInstance();
    
    Cimage*   poGetImage(const Cstring& strPath, bool bCollectGarbage = false);
    Cimage*   poGetImage(const int nDim0, const int nDim1, const eImage_data_types eD, Cstring& sName, bool bCollectGarbage = false);
    
    bool      bFreeImage(const Cstring& strPath, bool bDelete=true);
private:
    void vCollectGarbage();

    std::map<Cstring, Cimage*>    m_map_pImages;
    
    static CImagePool*  ms_poInstance;
};
#endif // DT_IMAGE_POOL_H
