//
// Copyright (c) 2005 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CoreSocket.cpp       Initial author: John Edwards  13-May-2005
// This file contains the definitions of class CDTCoreSocket

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

#ifndef CORESOCKET_H
#define CORESOCKET_H

#include "pinet.h"  // ptypes stuff
#include "Cstring.h"

typedef MSCString CString;

class DTREK_EXPORT CDTCoreSocket  
{
public:
     // public member vairables
     bool m_bConnected;

private:
     // private member variables
     const char * m_strHostIP;
     int m_iPort;
     bool m_bServer;
     CString m_strLastError;

     pt::ipstream * m_hConnectionSocket;
     pt::ipstmserver * m_hListenSocket;

public:
     int IsConnected();
     // public methods
     const char * GetLastError() const { return m_strLastError; }
     int ConnectToServer(const int iPort, const char * hostIP);
     bool SendString(CString str);
     bool ReceiveString(CString& strReceivedString);

    // These methods are private to enforce proper usage from within d*TREK.
private:
     bool AcceptClientConnection(const int iPort, const int TimeOut = 10);
//     CDTCoreSocket* AcceptMultipleClientConnections( const int iPort );
     bool Init(const char * hostIP, const int iPort, bool bServer);
     bool SetMaxPendingConnections(int iMaxConn);
     void Disconnect();

     bool ReceiveDataBlock(char * pData, const int iMax, int& BytesRead);
     int DataAvailable();
     bool SendDataBlock(const char * lpData, const int iDataSize);

public:
     // constructors and destructors
     CDTCoreSocket(const char * hostIP, const int port, bool bServer);
     CDTCoreSocket();
     virtual ~CDTCoreSocket();

protected:
     // protected methods
     pt::ipstmserver * CreateServerSocket(int iPort);

};

#endif // !CORESOCKET_H
