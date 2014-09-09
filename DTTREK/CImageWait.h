//
// Copyright (c) 2006 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CImageWait.h    Initial author: RB   07-Feb-2006
// This file contains the definitions of class CImageWait

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
#ifndef DT_IMAGE_WAIT_H
#define DT_IMAGE_WAIT_H

#include "dtreksys.h"

#include "Cstring.h"

class CImageWait;

class DTREK_EXPORT CImageWait
{
public:
    CImageWait();
    
    CImageWait(const CImageWait& oWait);
    
    ~CImageWait();

    CImageWait&  operator=(const CImageWait& oWait);

public:
    void vSetWaitTime(const double fTime){m_fWaitTime = fTime;}
    
    bool bWaitRequired()const{return m_fWaitTime > 0.0;}

    bool bDoWaitForFirstImage(Cstring& sImagePath);
    bool bDoWaitForNextImage(Cstring& sImagePath);
    bool bSaveImageSize(Cstring& sImagePath);

private:
    double      m_fWaitTime; // Max number of seconds to wait for an image
    
    int         m_nFileSize;
    int         m_nPercentFilesize;
};

#define DTREK_IMAGE_WAIT_GET        0x0001
#define DTREK_IMAGE_WAIT_DELETE     0x0002
CImageWait* g_poGetImageWaitObjectPtr(DTREK_WORD wCtrl);

#define THE_IMAGE_WAIT_PTR              (g_poGetImageWaitObjectPtr(DTREK_IMAGE_WAIT_GET)) 
#define DELETE_THE_IMAGE_WAIT_PTR       (g_poGetImageWaitObjectPtr(DTREK_IMAGE_WAIT_DELETE))

void g_vSetImageWaitTime(const double fTime);

#endif   // !DT_IMAGE_WAIT_H
