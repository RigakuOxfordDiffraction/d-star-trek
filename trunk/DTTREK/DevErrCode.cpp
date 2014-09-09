//
// Copyright 1997 Molecular Structure Corporation
//                9009 New Trails Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// DevErrCode.cc    Initial author: T.L.Hendrixson           23-Mar-1998
//  This file contains error code routines.
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

#include "DevErrCode.h"
#include "Cstring.h"

Cstring sGetDevErrString(int nCode)
{
   Cstring s;

   if(DEV_SUCCESS == nCode)
      s = "Success";
   else if(DEV_NOTCONNECTED == nCode)
      s = "Not connected, communications error";
   else if(DEV_INVALIDARG == nCode)
      s = "Invalid argument passed to routine";
   else if(DEV_INVALIDMODE == nCode)
      s = "In invalid mode";
   else if(DEV_INVALIDSTATE == nCode)
      s = "In invalid state";
   else if(DEV_INVALIDSETPOINT == nCode)
      s = "At invalid set point";
   else if(DEV_TIMEOUT == nCode)
      s = "Timeout";
   else if(DEV_IOCFAILED == nCode)
      s = "IOC hardware failure";
   else if(DEV_INVALIDAXIS == nCode)
      s = "Invalid axis specified";
   else if(DEV_FAILED == nCode)
      s = "General failure";
   else if(DEV_WRONGCOMMAND == nCode)
      s = "Wrong or unexpected command";
   else if(DEV_UNKNOWNERROR == nCode)
      s = "Unknown error";
   else if(DEV_INVALIDTYPE == nCode)
      s = "Invalid type specified";
   else if(DEV_COLLISION == nCode)
      s = "Collision will/has occurred";
   else if(DEV_WARNING == nCode)
      s = "Non-fatal warning";
   else if(DEV_ABORTED == nCode)
      s = "Action aborted";
   else if(DEV_INVALIDSYNTAX == nCode)
      s = "Invalid comand syntax";
   else
      s = "Unsupported error code";

   return s;
}

