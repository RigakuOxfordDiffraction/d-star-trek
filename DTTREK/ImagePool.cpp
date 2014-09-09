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
// ImagePool.cpp      Initial author: R.B.     06-June-2004
// This file is a storage of re-usable implementations for CDetector objects.
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
 
#include "ImagePool.h"

CImagePool* CImagePool::ms_poInstance = NULL;

CImagePool::CImagePool()
{
}

CImagePool::~CImagePool()
{
//    if( m_pInstance ) 
//    {
//        delete   m_pInstance;
//        m_pInstance = NULL;
//    }
}

CImagePool* CImagePool::poGetInstance()
{
    if( !ms_poInstance )
    {
        ms_poInstance = new CImagePool();
    }
    
    return ms_poInstance;
}

Cimage* CImagePool::poGetImage(const Cstring& strPath, bool bCollectGarbage)
{
    Cimage*     pImage = NULL;
    
    std::map<Cstring, Cimage*>::iterator         oIt;
    
    oIt = m_map_pImages.find(strPath);
    
    if( oIt != m_map_pImages.end() )
    {
        pImage = (*oIt).second;
    }
    else // create a new image
    {
        if( bCollectGarbage )
            vCollectGarbage();
        
        pImage = new Cimage(strPath);
        m_map_pImages.insert(std::pair<Cstring,Cimage*>(strPath,pImage));
    }

    pImage->nAddRef();
    return pImage;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cimage* CImagePool::poGetImage(const int nDim0, const int nDim1, const eImage_data_types eD, 
                               Cstring& strName, bool bCollectGarbage)
{
    Cimage*     pImage = NULL;
    
    std::map<Cstring, Cimage*>::iterator         oIt;
    
    oIt = m_map_pImages.find(strName);
    
    bool       bFound = false;
    if( oIt != m_map_pImages.end() )
    {
        pImage = (*oIt).second;
        
        if( pImage->nGetDimension(0) == nDim0 && pImage->nGetDimension(1) == nDim1 && pImage->nGetDataType() == (int)eD )
            bFound = true;
    }

    if( !bFound ) // create a new image
    {
        if( bCollectGarbage )
            vCollectGarbage();

        pImage = new Cimage(nDim0, nDim1, eD, strName);
        m_map_pImages.insert(std::pair<Cstring,Cimage*>(strName,pImage));
    }

    pImage->nAddRef();
    return pImage;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CImagePool::vCollectGarbage()
{
    // Clean the map of all images that are not in use
    Cimage*     pImage = NULL;
    
    std::map<Cstring, Cimage*>::iterator         oIt = m_map_pImages.begin();
    
    while( oIt != m_map_pImages.end() )
    {
        pImage = (*oIt).second;
        
        if( pImage && 0 == pImage->nGetRefCount() )
        {
            Cstring     strName = pImage->sGetName();
            delete pImage;
            pImage = NULL;
            m_map_pImages.erase(oIt++);
        }
        else 
            ++oIt;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Freeing means removing a reference and optionally deleting. If not deleting, make sure you collect "garbage", when
// creating a new image.
bool CImagePool::bFreeImage(const Cstring& strPath, bool bDelete)
{
    Cimage*     pImage = NULL;
    
    std::map<Cstring, Cimage*>::iterator      oIt;
    
    oIt = m_map_pImages.find(strPath);
    
    if( oIt != m_map_pImages.end() )
    {
        pImage = (*oIt).second;
        
        if( pImage )
        {
            int     nRef = pImage->nReleaseRef();
            
            if( bDelete && 0 == nRef )
            {
                delete  pImage;
                pImage = NULL;
                m_map_pImages.erase(strPath);
            }
        }
        else
            return false;  // image path (keyword) not found
    }
    else
        return false; // image path (keyword) not found
    
    return true;
}

