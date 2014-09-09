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
// DTMainCell.h     Initial author:                RB 10-Sep-2004

// Wrapper-class around former dtcell entry point

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

#ifndef DT_MAIN_CELL_H
#define DT_MAIN_CELL_H

#include "DTMain.h"

class DTREK_EXPORT CDTMainCell : public CDTMain
{
public:
    CDTMainCell();
public:
    virtual int nExecute(unsigned int argc, char* argv[]);
private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage);

    void vNotifyClientAboutFileWritten(Cstring& sText, DTREK_WORD wCtrl);

private:
                                    // Thad's comments:
    bool    m_bCheckCentric;        // Are we checking for centricity.  If not, then m_bIsCentric will be correctly set.  
    bool    m_bCheckLaue;           // Are we checking the Laue group.  If not, Claue::m_eLaue will be correctly set.
    bool    m_bSpacegroupCheck;     // Always test the spacegroup unless otherwise specified.
    bool    m_bReduce;              // Call cell reduction code.
};
#endif   // !DT_MAIN_CELL_H
