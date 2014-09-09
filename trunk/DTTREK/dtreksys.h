#ifndef DT_DTREKSYS_H
#define DT_DTREKSYS_H
//
// Copyright (c) 1996 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtreksys.h        Initial author: J.W. Pflugrath           17-Nov-1996
//    This file is the header file for some modules
//    which implement some system dependent functions
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

#ifdef WIN32
    #ifdef _DEBUG
        #define _CRTDBG_MAP_ALLOC
        #include <stdlib.h>
        //#include <crtdbg.h>
    #endif
#endif

#include "Dtrek.h"
#include "dtrekdefs.h"
#include "Cstring.h"
#include "dtarray.h"


#include <algorithm> // for std::swap

//+Definitions and constants

enum eCPU_byte_orders {
  eCPU_other_endian,
  eCPU_big_endian,
  eCPU_little_endian
};

// Function prototypes

#ifdef WIN32

DTREK_EXPORT int nint(double d);
//DTREK_EXPORT long rint(double d);

#endif

//Cstring GetInstallDirectory(void);
DTREK_EXPORT Cstring sTransSymbol(const Cstring& rsStringIn);
DTREK_EXPORT Cstring sGetEnv(const Cstring& rsStringIn);
DTREK_EXPORT bool bGetEnv(const Cstring& rsStringIn, double& dValue);

DTREK_EXPORT int   nPutEnv(const Cstring& rsKeyword, const Cstring& rsValue);
DTREK_EXPORT int   nSetCWD(const Cstring& rsNewCWD);
DTREK_EXPORT Cstring sGetCWD(void);
DTREK_EXPORT Cstring sGetTime(void);
DTREK_EXPORT Cstring sGetDate(void);

DTREK_EXPORT double dGetTimeOfDay(void);
DTREK_EXPORT int nWait(int nMillSecond);

DTREK_EXPORT int  nDoSystemCommand(const Cstring& rsCommand, int *pnPID=NULL);
DTREK_EXPORT int  nDoSystemCommand(const Cstring& rsCommand, FILE **ppSubprocessStdIn, int *pnPID=NULL);
#ifdef SSI_PC
DTREK_EXPORT int  nDoSystemCommand2(const Cstring& rsCommand, const Cstring& rsOutFile);
#endif
DTREK_EXPORT int  nFileRename(const Cstring& rsFilenameOld, const Cstring& rsFilenameNew);
DTREK_EXPORT int  nFileCopy(const Cstring& rsFilenameOld, const Cstring& rsFilenameNew);
DTREK_EXPORT bool bFileExists(const Cstring& rsFilename);
DTREK_EXPORT int  nFileAppendVersion(const Cstring& rsFilename, const bool bCheckEnv);
DTREK_EXPORT int  nFileCreateDirectory(const Cstring& rsDirName);
DTREK_EXPORT int  nFileDelete(const Cstring& rsFilename);

DTREK_EXPORT double dGetFreeSpaceOnDisk(Cstring &sFileOnDisk);
DTREK_EXPORT double dGetFreeSpaceOnDisk(const char *pcFileOnDisk = "");

DTREK_EXPORT Cstring sFileBuildName(Cstring &rsDirectory, Cstring &rsFilename);
DTREK_EXPORT Cstring sFileBuildName(char *pcDirectory, Cstring &rsFilename);
DTREK_EXPORT Cstring sFileBuildName(Cstring &rsDirectory, char *pcFilename);
DTREK_EXPORT Cstring sFileBuildName(char *pcDirectory, char *pcFilename);
DTREK_EXPORT Cstring sFileGetDirectory(const Cstring& rsFilename);
DTREK_EXPORT Cstring sFileGetBasename(const Cstring &rsFilename);
DTREK_EXPORT Cstring sFileGetExtension(const Cstring &rsFilename);
DTREK_EXPORT Cstring sFileGetTempName(const char *dir, const char *pfx);
DTREK_EXPORT Cstring sFileGetFullPath(Cstring& sRelPath);
DTREK_EXPORT Cstring sDtrekGetPrefix(void);
DTREK_EXPORT Cstring sBuildScanTemplate(const Cstring &rsFilename, 
				 const int nPlaces=4, int *pnSeqNum=NULL);

DTREK_EXPORT Cstring sBuildStrategyScanPrefix(const Cstring& sStrategyPrefix, const Cstring& sPre, int nScanNumber);
DTREK_EXPORT Cstring sBuildScanPrefix(const Cstring& sPre, int nScanNumber);

DTREK_EXPORT Cstring sDtrekGetModuleName(void);
DTREK_EXPORT void  vDtrekSetModuleName(const Cstring& rsString);
DTREK_EXPORT void  vPrintCopyrightInfo(void);

DTREK_EXPORT int    nGetByteOrderCPU(void);
DTREK_EXPORT Cstring sGetByteOrderCPU(void);

DTREK_EXPORT LONG  lFileGetSize(const Cstring& rsFilename);
DTREK_EXPORT int nKillProcess(LONG lPID, int nSignal);

DTREK_EXPORT Cstring sGetCommandLine(int argc, char *argv[], const int nMaxWidth=72);
DTREK_EXPORT int nModifyCommandLine(unsigned int& nArgc,char*** pargv);
DTREK_EXPORT Cstring sTabFormat(Cstring& sString, int nTabCount,int nMaxWidth = 72);

DTREK_EXPORT int nGetBatchNumber(Cstring& sBatch,Cstring& sScanID);
DTREK_EXPORT void vMakeBatchName(int nBatchNumber, const char* pcScanPrefix, Cstring& sBatchName);
DTREK_EXPORT void vModifyDISPLAY(void);
DTREK_EXPORT int nTokenize(Cstring& sStringIn, 
                           std::vector<Cstring>& asTokensOut,
                           std::vector<Cstring>* pasSymbols = NULL);

const double g_fTableAverage = 1.123456e-90;
const double g_fTableMax = 2.123456e-90;
const double g_fTableMin = 3.123456e-90;
const double g_fTableSum = 4.123456e-90;

DTREK_EXPORT int nPrintTable(Cstring& sTableTitle,
                             std::vector<Cstring>& asTableHeaders,
                             std::vector<Cstring>& TableFormat,
                             itr<double>& afValues,
                             bool bLastRowIsSummary);

extern Cstring* g_psPrintfLog;
DTREK_EXPORT int lprintf(const char* pcFormat,...);

#ifdef VC9
char* dt_strcpy(char *strDestination, const char *strSource);
char* dt_strncpy(char *strDest, const char *strSource, size_t count);
char* dt_strtok(char *strToken, const char *strDelimit);
#endif

extern Cstring g_sCommandLine;

float fSwapFloat(float fFloat);
float fVAXtoIEEE(float *pfTemp);
float fIEEEtoVAX(float *pfTemp);

inline INT4
nSwapLong(UINT4 ulLong)
{
  return ( (ulLong << 24) | ((ulLong << 8) & 0x00ff0000) |
	          ((ulLong >> 8) & 0x0000ff00) | (ulLong >> 24));
}

#define DTREK_MAX_LINE_LENGTH        256
void set_debug_file_path(Cstring& sPath);
void set_debug_env_var_name(Cstring& sVar);
void print_debug_string(Cstring& s);
void print_debug_string(char* pc);
void create_debug_file();
void close_debug_file();
void delete_debug_file();
struct ifreq* DTREK_next_ifreq(struct ifreq* req);
DTREK_EXPORT int  nDTREKGetIPAddress(Cstring *psIPAddress, int a4nIPField[4]);

#define REPORT_ERROR_FILE_LINE {  printf("Bad Logic Line %d of file '%s'\n", __LINE__ , __FILE__ ); GETCH; exit(0); };

#endif   // DT_DTREKSYS_H
