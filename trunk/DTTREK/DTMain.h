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
// DTMain.h     Initial author:  RB 15-Mar-2004

// Base wrapper class around former dt main entry points

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

#ifndef DT_MAIN_H
#define DT_MAIN_H

#define DTREK_ERROR(code,msg) {vError(code,msg);return (code);}

#ifdef WIN32 
#ifndef SSI_PC
//#include "Forcelib.h"
#endif
#endif

#include "Dtrek.h"
#include "dtreksys.h"

#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#define getline(x,y) y = CCrclHelper::GetInstance()->SendTclCmd(x)
#endif

#include "Cimage_header.h"

#include "DTCoreSocket.h"


class DTREK_EXPORT CDTMain 
{
public:
    CDTMain();
    CDTMain(unsigned int argc, char* argv[]);

    virtual ~CDTMain();

protected:
    bool    m_bExternalHeader;

private:
    void vInitialize();
    void vClearCommandArguments();

public:
#define     DTMAIN_ERROR_UNRECOGNIZED_COMMAND_OPTION    7
    virtual int nExecute(unsigned int argc, char* argv[])=0;
    
    //virtual int nExecute()=0; // a stub for the future... nExecute shouldn't be receiving argc and argv as arguments... 

    void vSetExternalHeader(Cimage_header* pHeader, bool bDeleteCurrentHeader=true);

protected:
    Cimage_header*      m_poHeader;
    
    Cstring             m_sHeaderOut;
    Cstring             m_sReflnListOut;
    
    Cstring             m_strCommandLineOptions;
    int                 m_nDTREKCommandArgumentCount;
    char**              m_ppcDTREKCommandArguments;

	// A socket object used for communicating status back to
    // clients.  This is NOT intended for input into d*TREK.
    // We may want to use a wrapper class to enforce this.
    CDTCoreSocket*        m_poSocket;

    bool bParseCalculateResolution(Cstring& sResoCommandLineOption, float& fResoMin, float& fResoMax, Cstring& strError);
    bool bParseCalculateResolution(Cstring& sResoCommandLineOption, double& dResoMin, double& dResoMax, Cstring& strError);

    bool bCreateHeader(Cstring& sFileName);
    void vCreateHeader();
    void vCommandLineToCommandArgumentArray(char* strCommand);
    void vCommandLineToCommandArgumentArray(unsigned int argc, char* argv[]);
    void vGetCommandArrayString(Cstring& sCommand);


    bool bCalculateDetectorResolution(float& fFirstResolution, 
                                      float& fSecondResolution, 
                                      float& fResolutionMaxEdge,
                                      Cstring& sError);
    
    bool bCalculateDetectorResolution(double& dFirstResolution, 
                                      double& dSecondResolution, 
                                      double& dResolutionMaxEdge,
                                      Cstring& sError);

    void vSetCommandArguments(unsigned int argc, char* argv[]);
    bool bIsHelpRequest();
    
    bool bGetCommandOption(const char* ccOption, Cstring& sArgs, bool bRemove=true, int* piStart=NULL, int* piEnd=NULL);
    bool bGetCommandOption(const char* ccOption, bool bRemove=true, int* piStart=NULL, int* piEnd=NULL);
    
    bool bGetCommandOption(const char* ccOption, 
                           std::vector<double>& daArgs, 
                           char* pcListSeparators, 
                           bool bRemove=true, int* piStart=NULL, int* piEnd=NULL);

    bool bGetCommandOption(const char* ccOption, 
                           std::vector<Cstring>& asArgs, 
                           char* pcListSeparators, 
                           bool bRemove=true, int* piStart=NULL, int* piEnd=NULL);
    
    bool bGetCommandOption(const char* ccOption, double& dArg, bool bRemove=true, int* piStart=NULL, int* piEnd=NULL);

    bool bRemoveCommandOption(const char* ccOption);
    bool bRemoveCommandArgs(const int iStart, const int iEnd);
    bool bUpdateCommandOption(const char* ccOption, Cstring& sArgs);
    void vPrependCommandOption(const char* ccOption);

    bool bOpenSocket(int nPort);

    void vSetOutputHeaderPath();
    void vSetOutputReflnListPath();

private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage)=0;
};
#endif   //!DT_MAIN_H

