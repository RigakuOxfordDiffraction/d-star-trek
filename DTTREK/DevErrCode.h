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
// DevErrCode.h     Initial author: Bob Daly
//
// The original version of this file was developed for the CDevxxx d*TREK
// classes which involve communication with an EPICS IOC.
// This file has been expanded to contain additional error codes and also
// a routine returning an string representation of an error code.
//
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

#ifndef DEV_ERROR_CODE_H
#define DEV_ERROR_CODE_H

#include "Cstring.h"

// success
#define DEV_SUCCESS         0  /* dev success  				*/
// errors
#define DEV_NOTCONNECTED    1  /* communications error to IOC		*/
#define DEV_INVALIDARG      2  /* invalid argument passed to dev calls	*/
#define DEV_INVALIDMODE     3  /* in wrong mode after dev call		*/
#define DEV_INVALIDSTATE    4  /* in wrong mode after dev call		*/
#define DEV_INVALIDSETPOINT 5  /* in wrong mode after dev call		*/
#define DEV_TIMEOUT         6  /* time out 				*/
#define DEV_IOCFAILED       7  /* IOC hardware failed			*/
#define DEV_INVALIDAXIS     8  /* invalid axis specified */
#define DEV_FAILED          9  /* general failure */
#define DEV_WRONGCOMMAND   10  /* wrong command received */
#define DEV_UNKNOWNERROR   11  /* unknown or unexpected error */
#define DEV_INVALIDTYPE    12  /* invalid hardware type specified */
#define DEV_COLLISION      13  /* in/through gonio collision area */
#define DEV_WARNING        14  /* non-fatal error/warning */
#define DEV_ABORTED        15
#define DEV_EVENT          16  /* event callback signaled return */
#define DEV_INVALIDSYNTAX  17  /* invalid command syntax */

Cstring sGetDevErrString(int nCode);

#endif  /* DEV_ERROR_CODE_H */
