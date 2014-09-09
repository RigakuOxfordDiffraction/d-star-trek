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
// This file contains the member functions of class CDTCoreSocket

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

#include <stdio.h>
#include "DTCoreSocket.h"

using namespace pt;

// TODO : Consider removing m_bConnected state variable if it is
//        redundant since m_hConnectionSocket == null also represents not connected

static cset VALID_CHARS = cset("*") - cset("~00");

CDTCoreSocket::CDTCoreSocket()
{
     // initialize state variables
     m_bServer = false;
     m_bConnected = false;
     m_iPort = 0;
     m_strHostIP = "127.0.0.1";

     // init socket handles
     m_hConnectionSocket = NULL;
     m_hListenSocket = NULL;
}

CDTCoreSocket::CDTCoreSocket(const char * hostIP, const int iPort, bool bServer)
{
     // initialize state variables
     m_bServer = bServer;
     m_bConnected = false;
     m_iPort = iPort;
     m_strHostIP = hostIP;

     // init socket handles
     m_hConnectionSocket = NULL;
     m_hListenSocket = NULL;
}

CDTCoreSocket::~CDTCoreSocket()
{
     // close open socket handles
     if (m_hConnectionSocket != NULL)
     {
          m_hConnectionSocket->close();
	  delete m_hConnectionSocket;
     }
     if (m_hListenSocket != NULL)
     {
	  delete m_hListenSocket;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   SendString
// 
// Parameters
//     lpString     -  String to be sent over the stream connection
// 
// Returns
//    true if successful
//     false if failure.
//
//
// Preconditions
//   valid connection has been established
//        
// Postconditions
//   data has been sent
//        
//*************************************************************************************
bool CDTCoreSocket::SendString(CString str)
{
     if (!IsConnected())
     {
          m_strLastError = "Socket not connected";
          return false;
     }

     // write data to socket
     try
     {
          m_hConnectionSocket->put(str);
          m_hConnectionSocket->put('\0');
          m_hConnectionSocket->flush();
          return true;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   SendSendDataBlock
// 
// Parameters
//     lpData  - Data to be sent over the connection
//     iDataSize    - size of the data to be sent
// 
// Returns
//    true if successful
//     any other if failure.
//
//
// Preconditions
//   Valid connection has been established
//        
// Postconditions
//   data has been sent
//        
//*************************************************************************************
bool CDTCoreSocket::SendDataBlock(const char * lpData, const int iDataSize)
{
     if ( ( !lpData ) || ( !iDataSize ) )
          return false;

     if (!IsConnected())
     {
          m_strLastError = "Socket not connected";
          return false;
     }

     // write data to socket
     try
     {
          m_hConnectionSocket->write(lpData, iDataSize);
          m_hConnectionSocket->flush();
          return true;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   ReceiveString
// 
// Parameters
//        strReceivedString - reference to a CString that will receive the incomming command
// 
// Returns
//    true if successful
//     false if failure.
//
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
bool CDTCoreSocket::ReceiveString(CString& strReceivedString)
{
     if (!IsConnected())
     {
          m_strLastError = "Socket not connected";
          return false;
     }
     
     try
     {
          strReceivedString = m_hConnectionSocket->token(VALID_CHARS);
          m_hConnectionSocket->get(); // skips the null terminator
          return true;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   DataAvailable
// 
// Parameters
//     None
// 
// Returns
//    -1 = not connected
//      0 = Connected, no input
//      1 = Connected, input pending (data available)
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
int CDTCoreSocket::DataAvailable()
{
     if (!IsConnected())
          return -1;
     if (m_hConnectionSocket->waitfor(0)) 
          return 1;
     else
          return 0;
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   ReceiveDataBlock
// 
// Parameters
//        pData          - buffer for incomming data
//        iDataSize - sized of pData buffer
//        BytesRead - number of bytes placed into buffer during read
// 
// Returns
//    1 if successful
//     any other if failure.
//
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
bool CDTCoreSocket::ReceiveDataBlock(char * pData, const int iMax, int& BytesRead)
{
     if (!IsConnected())
     {
          m_strLastError = "Socket not connected";
          return false;
     }
     
     try
     {
          int count = m_hConnectionSocket->read(pData, iMax);
          BytesRead = count;
          return true;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   SetMaxPendingConnections
// 
// Parameters
//        iMaxConn  - maximum number of incomming connection request that will be queued
// 
// Returns
//    1 if successful
//     any other if failure.
//
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
bool CDTCoreSocket::SetMaxPendingConnections(int iMaxConn)
{
     // This function no longer does anything useful but has been maintained for
     // backwards compatibility.  The truth is that in simple client/server architectures
     // controlling the max number of pending connections serves no purpose.  Also
     // even if this were to default to only 1, it would have any affect on a server
     // that wants to accept multiple client connections.  The intended use of
     // multiple pending connections is in a high traffic network where multiple 
     // clients are competing for a socket connection, which is not the case here.

     if (m_hListenSocket == NULL)
          return false;
     return true;
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   AcceptClientConnection
// 
// Parameters
//        iPort     - port to listen for connections on
// 
// Returns
//    true if successful
//     false if failure.
//
//
// Preconditions
//  valid connection must not exist
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
bool CDTCoreSocket::AcceptClientConnection(const int iPort, const int TimeOut)
{
     if (IsConnected())
     {
          m_strLastError = "Socket already connected";
          return false;
     }

     // create socket to use to listen for incomming connections
     if (m_hListenSocket == NULL)
     {
          m_iPort = iPort;
          m_bServer = true;

          m_hListenSocket = CreateServerSocket(iPort);
     }

     // accept a client connection
     ipstream * client = NULL;
     try
     {
          client = new ipstream();
          bool success = m_hListenSocket->serve(*client, -1, TimeOut);

          if (success) 
          {
               m_hConnectionSocket = client;
               m_bConnected = true;
          }
          return success;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          delete client;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   ConnectToServer
// 
// Parameters
//        iPort     - server port to connect to
//        hostIP    - IP address or DNS name of server to connect to
// 
// Returns
//    true if successful
//     false if failure.
//
//
// Preconditions
//   socket has been initialized
//        
// Postconditions
//   connection to a server established
//        
//*************************************************************************************
int CDTCoreSocket::ConnectToServer(const int iPort, const char * hostIP)
{
     m_bServer = false;
     if (hostIP)
     {
          m_strHostIP = hostIP;
     }
     m_iPort = iPort;

     // resolve IP address
     ipaddress addr = phostbyname(hostIP);
     if (addr == ipnone) {
          m_strLastError = "Host not found [";
          m_strLastError += hostIP;
          m_strLastError += "]";
          return false;
     }

     try 
     {
          // connect to server
          ipstream * connection = new ipstream(addr, iPort);
          connection->open();
          // save connection info
          m_hConnectionSocket = connection;
          m_bConnected = true;
          return true;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          return false;
     }
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   Disconnect
// 
// Parameters
//     None
// 
// Returns
//    1 if successful
//     any other if failure.
//
//
// Preconditions
//   Valid connection exists
//        
// Postconditions
//   Connected closed and socket handle freed.
//        
//*************************************************************************************
void CDTCoreSocket::Disconnect()
{
     // if connected then initiate shutdown sequence
     if (m_bConnected == true )
     {
          if (m_hConnectionSocket != NULL)
          {
               m_hConnectionSocket->close();
               delete m_hConnectionSocket;
               m_hConnectionSocket = NULL;
          }
          m_bConnected = false;
      }
}


//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   Init
// 
// Parameters
//        hostIP    - ip address or DNS name of server
//        iPort     - server port number
//        bServer - server or client
// 
// Returns
//    1 if successful
//     any other if failure.
//
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
bool CDTCoreSocket::Init(const char * hostIP, const int iPort, bool bServer)
{
     m_bServer = bServer;
     m_strHostIP = hostIP;
     m_iPort = iPort;

     if ( bServer == true )
     {
          // Create socket to listen for connections
          m_hListenSocket = CreateServerSocket( m_iPort );
     }

     return true;
}

//*************************************************************************************
// Class:      CDTCoreSocket
// Function:   CreateServerSocket
// 
// Parameters
//        iPort     - port number for server to listen on
// 
// Returns
//    1 if successful
//     any other if failure.
//
//
// Preconditions
//   No other members of this object have been called
//        
// Postconditions
//   Initialization completed successfully
//        
//*************************************************************************************
ipstmserver * CDTCoreSocket::CreateServerSocket(int iPort)
{
     ipstmserver * server;
     try
     {
          server = new ipstmserver();               
          server->bindall(iPort);
          return server;
     }
     catch (estream * e)
     {
          m_strLastError = e->get_message();
          delete e;
          delete server;
          return NULL;
     }
}

// test if socket is connected
int CDTCoreSocket::IsConnected()
{
     if (m_hConnectionSocket == NULL) 
          return false;

     int status = m_hConnectionSocket->get_status();
     if (status == IO_CREATED || status == IO_EOF ||
         status == IO_CLOSING || status == IO_CLOSED) 
     {
          return false;
     }

     return true;
}
