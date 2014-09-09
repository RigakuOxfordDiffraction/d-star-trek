//
// Copyright (c) 1998 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DevEnumeration.h      Initial author: tlh    JWP           Jul 1998
//  This file contains the various enumerations for the ADev classes, 
//     as well as character string representations of those enumerations.
//     All of the enumerations were placed in this file to give non-derived
//     classes access to the enumerated values without having to link
//     to the compiled classes (which one was forced to do since the 
//     ADev classes were (naturally) defined in the ADev*.h files).
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
//+Description
//
//+ToDo
//
#ifndef DEV_ENUMERATION_H
#define DEV_ENUMERATION_H

// Enumerations for various states and modes

// DevCounter state

enum eDevCounter_State
{
    eDevCounterState_NotAvailable,
    eDevCounterState_Unknown,
    eDevCounterState_Initializing,
    eDevCounterState_Counting,
    eDevCounterState_Idle,
    eDevCounterState_NumStates
};

// DevDetector mode

enum eDevDetector_Mode
{
    eDevDetMode_StillClosed,        // These modes might be better suited to
    eDevDetMode_StillOpen,          // the DevScan class...
    eDevDetMode_ScanClosed,
    eDevDetMode_ScanOpen,
    eDevDetMode_NumModes
};

// DevDetector state

enum eDevDetector_State
{
    eDevDetState_Unknown,
    eDevDetState_NotReady,
    eDevDetState_Ready,
    eDevDetState_Exposing,
    eDevDetState_Reading,
    eDevDetState_Storing,
    eDevDetState_Erasing,
    eDevDetState_Done,
    eDevDetState_Initializing,
    eDevDetState_Initialized,
    eDevDetState_Moving,
    eDevDetState_Failed,
    eDevDetState_FailedRetry,
    eDevDetState_NumStates
 };

// DevGoniometer mode

enum eDevGoniom_Mode
{
    eDevGoniomMode_Unknown,
    eDevGoniomMode_Slew,
    eDevGoniomMode_Scan,
    eDevGoniomMode_NumModes
};

// DevGoniometer state

enum eDevGoniom_State
{
    eDevGoniomState_Unknown,
    eDevGoniomState_Initializing,
    eDevGoniomState_NotReady,
    eDevGoniomState_Ready,
    eDevGoniomState_Moving,
    eDevGoniomState_Done,
    eDevGoniomState_Aborted,
    eDevGoniomState_TimedOut,
    eDevGoniomState_Failed,
    eDevGoniomState_Collided,
    eDevGoniomState_Jammed,
    eDevGoniomState_Initialized,
    eDevGoniomState_NumStates
};

// DevScan mode

enum eDevScan_Mode
{
    eDevScanMode_Unknown,                         //For multiple DARK images
    eDevScanMode_StillClosed,                     //For multiple DARK images
    eDevScanMode_StillOpen,                       //For Multiple STILLS
    eDevScanMode_ScanClosed,                      //For Testing Purposes only
    eDevScanMode_ScanOpen,                        //For X-ray Diffraction Data
    eDevScanMode_MINOR,
    eDevScanMode_INVALID,
    eDevScanMode_NumModes
};

// DevScan state

enum eDevScan_State
{
    eDevScanState_Unknown,
    eDevScanState_NotReady,
    eDevScanState_Ready,
    eDevScanState_Starting,
    eDevScanState_InProgress,
    eDevScanState_Reading,
    eDevScanState_Pausing,
    eDevScanState_Paused,
    eDevScanState_Resuming,
    eDevScanState_Done,
    eDevScanState_Aborted,
    eDevScanState_NotTimedOut,
    eDevScanState_Failed,
    eDevScanState_Erasing,
    eDevScanState_Initializing,
    eDevScanState_NumStates
};

// DevShutter state

enum eDevShutter_State
{
    eDevShutterState_Unknown,
    eDevShutterState_Initializing,
    eDevShutterState_Open,
    eDevShutterState_Closed,
    eDevShutterState_Opening,
    eDevShutterState_Closing,
    eDevShutterState_TimedOut,
    eDevShutterState_Failed,
    eDevShutterState_NumStates
};

// DevSource state

enum eDevSource_State
{
    eDevSourceState_Unknown,
    eDevSourceState_Initializing,
    eDevSourceState_Available,
    eDevSourceState_NotAvailable,
    eDevSourceState_InProgress,
    eDevSourceState_NumStates
};


// Character string representations of the states and modes

static const char *g_pcCounterStates[] = 
{
    "Not Available",
    "Unknown",
    "Initializing",
    "Counting",
    "Idle",
    "Number of States"
};

static const char *g_pcDetectorModes[] =
{
    "Still Closed",
    "Still Open",
    "Scan Closed",
    "Scan Open",
    "Number of Modes"
};

static const char *g_pcDetectorStates[] = 
{
    "Unknown",
    "Not Ready",
    "Ready",
    "Exposing",
    "Reading",
    "Storing",
    "Erasing",
    "Done",
    "Initializing",
    "Initialized",
    "Moving",
    "Failed",
    "Failed Retry",
    "Number of States"
};

static const char *g_pcGoniomModes[] =
{    "Unknown",
     "Slew",
     "Scan",
     "Number of Modes"
};

static const char *g_pcGoniomStates[] =
{    "Unknown",
     "Initializing",
     "Not Ready",
     "Ready",
     "Moving",
     "Done",
     "Aborted",
     "Timed Out",
     "Failed",
     "Collided",
     "Jammed",
     "Initialized",
     "Number of States"
};

static const char *g_pcScanModes[] = 
{
    "Unknown",
    "Still Closed",
    "Still Open",
    "Scan Closed",
    "Scan Open",
    "MINOR",
    "INVALID",
    "Number of Modes"
};

static const char *g_pcScanStates[] = 
{
    "Unknown",
    "Not Ready",
    "Ready",
    "Starting",
    "In Progress",
    "Reading",
    "Pausing",
    "Paused",
    "Resuming",
    "Done",
    "Aborted",
    "Not Timed Out",
    "Failed",
    "Number of States"
};

static const char *g_pcShutterStates[] =
{    "Unknown",
     "Initializing",
     "Open",
     "Closed",
     "Opening",
     "Closing",
     "Timed Out",
     "Failed",
     "Number of States"
};

static const char *g_pcSourceStates[] =
{
    "Unknown",
    "Initializing",
    "Available",
    "Not Available",
    "In Progress",
    "Number of States"
};

#endif /* DEV_ENUMERATION_H */
