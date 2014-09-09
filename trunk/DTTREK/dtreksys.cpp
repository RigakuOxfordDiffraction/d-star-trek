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
// dtreksys.cc     Initial author: J.W. Pflugrath           17-Nov-1996
//    This file contains some system dependent d*TREK functions
//    that do not belong to any class.
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

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "dtarray.h"

// May need to #ifdef for various operating systems

#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>

#ifdef WIN32
#include <winsock2.h>
#include <windows.h>           // for defs of things like DWORD
#include <direct.h>
#include <sys/timeb.h>
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/param.h>
#ifdef OSF1
#include <sys/mount.h>
#elif !defined(OSX)
#include <sys/statfs.h>
#endif
#endif

#ifdef LINUX
#include <sys/timeb.h>
#include <limits.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <net/if.h>

#endif

#ifdef OSX
#include <sys/param.h>
#include <sys/mount.h>
#endif

#ifdef SUNOS
// This gets the SunOS swap routine
#include <algorithm>
#endif

#include "stdarg.h"
#include "dtreksys.h"           // Function prototypes


#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
#endif

//+Definitions, constants, and initialization of static member variables

const Cstring Gs_DTREK_OVERWRITE = "DTREK_OVERWRITE";

//+Code begin

//+Public functions

Cstring sTransSymbol(const Cstring& rsStringIn)
{
  // Translate any embedded symbols (aka environment variables) in a string.
  // In the input symbols take the form $(symbol).
  // They may be embedded anywhere in the input string.
  //
  // Suppose
  //    setenv HOME  home
  //    setenv TERM  term
  //    setenv HOMEterm nested
  // They may be nested:  $(HOME$(TERM))    ->   nested
  //                      $(TERM$(HOME))    ->   TERMhome
  //      or sequential:  $(HOME)$(TERM)    ->   hometerm
  //          or absent:  HOMETERM          ->   HOMETERM
  //            or null:  $()               ->
  //           or wrong:  $(HOME            ->   $(HOME
  //           or wrong:   HOME)            ->   HOME)
  // Unlike the unix getenv which returns NULL if a symbol
  // does not exist, this routine will return the symbol name,
  // i.e. the input string without the delimiters $( and ).

  int     nIndexL, nIndexR;
  static  Cstring s_sDollarLeft = "$(";
  static  Cstring s_sRight      = ")";
  Cstring sResult;
  Cstring sTemp;
  Cstring sTemp1;
  Cstring sRight;

  nIndexR = -1;
  sResult = rsStringIn;
  nIndexL = sResult.find(s_sDollarLeft);
  if (0 <= nIndexL)
    {
      // Potential environment variable found, so ...
      // Look for a nested environment variable with a little recursion
        
      sRight  = sTransSymbol(sResult.from(nIndexL+2));
        
      nIndexR = sRight.find(s_sRight);
      if (0 <= nIndexR)
        {
          // Found ), so translate what is in between
        
          sTemp  = sRight.before(nIndexR);
          sTemp1 = sGetEnv(sTemp);
          if ("" != sTemp1)
            {
              sTemp = sTemp1;
            }
          sResult = sResult.before(nIndexL) + sTemp + sRight.after(nIndexR);
        }
      else
        {
          cout << "WARNING in sTransSymbol: no ')' found!\n";
        }
    }
  return (sResult);
}

Cstring sGetEnv(const Cstring &rsStringIn)
{
  char *pcGetenv;
  pcGetenv = getenv(rsStringIn.string());   // Very UNIX!
  if (NULL == pcGetenv)
    return ("");
  else
    return ((Cstring)pcGetenv);
}
/////////////////////////////////////////////////////////////
bool bGetEnv(const Cstring& rsStringIn, double& dValue)
{
    Cstring     strValue = sGetEnv(rsStringIn);

    if( "" == strValue )
        return false;

    if( !strValue.bIsNumeric() )
        return false;

    dValue = atof(strValue.string());

    return true;
}
/////////////////////////////////////////////////////////////
int nPutEnv(const Cstring &rsKeyword, const Cstring& rsValue)
{
  // Add to or modify the environment

  Cstring *psTemp;

  // The following line introduces a memory leak, but putenv requires a
  // static variable, so this is an easy way to do it.
  // TODO: use sTransSymbol

  if ("" != rsValue)
    psTemp = new Cstring(rsKeyword + "=" + rsValue);
  else
    psTemp = new Cstring(rsKeyword + "=" + """");

#ifdef VC9
  return(_putenv(psTemp->string()));
#else
  return(putenv(psTemp->string()));
#endif

}

Cstring sGetCWD(void)
{
  // Return the value of the current working directory ready for prepending
  // to a filename. (May have slash, backslash or right bracket as last char!)

  Cstring sResult;
  char *pcCWD = NULL;
  int  nSize  = 256;
  int nLen;

#ifdef WIN32
  pcCWD = _getcwd((char *)NULL, nSize);
#else
  pcCWD = getcwd((char *)NULL, nSize);
#endif

  if (NULL != pcCWD)
    {
      sResult = pcCWD;
#ifdef WIN32
      nLen = sResult.length();
      if('\\' != sResult.GetAt(sResult.length()-1)) // not root directory
         sResult += '\\';
#else
      sResult += '/';
#endif
      free(pcCWD);
    }
  else
      sResult = "";

    return (sResult);
}

int nSetCWD(const Cstring& rsNewCWD)
{
  Cstring sTemp;

  // Remove preceding and trailing whitespace from sTemp

  sTemp = rsNewCWD;
  while (  (0 < sTemp.length())
         &&
         (  (' '  == sTemp.GetAt(0))
          || ('\t' == sTemp.GetAt(0))) )
    sTemp = sTemp.after(0);
  while (  (0 < sTemp.length())
         &&
         (   (' '  == sTemp.GetAt(sTemp.length()))
          || ('\t' == sTemp.GetAt(sTemp.length())) ) )
    sTemp = sTemp.before((int)sTemp.length());

  if (0 == sTemp.length())
    {
      sTemp = "$(HOME)";
    }
  sTemp = sTransSymbol(sTemp);

#ifdef WIN32
  return (_chdir(sTemp.string()));
#else
  return (chdir(sTemp.string()));
#endif
}

Cstring sGetTime(void)
{
  // Return the current time of day in format HH:MM:SS

  Cstring sTime;
  char   *pcTime = NULL;
  time_t  tTime;

  time(&tTime);
  pcTime = ctime(&tTime);

  // 012345678 012345678 012345
  // Fri Sep 13 00:00:00 1986\n\0

  sTime = "??:??:??";
  if (NULL != pcTime)
    {
      sTime.SetAt(0, pcTime[11]);
      sTime.SetAt(1, pcTime[12]);
      sTime.SetAt(3, pcTime[14]);
      sTime.SetAt(4, pcTime[15]);
      sTime.SetAt(6, pcTime[17]);
      sTime.SetAt(7, pcTime[18]);
    }
  return (sTime);
}

Cstring sGetDate(void)
{
  // Return the current date in format DD-MMM-YYYY

  Cstring sDate;
  char   *pcTime = NULL;
  time_t  tTime;

  time(&tTime);
  pcTime = ctime(&tTime);  // Note pcTime need not be freed

  // 012345678 012345678 012345
  // Fri Sep 13 00:00:00 1986\n\0

  // I am building this string this way because the Microsoft compiler
  // treats "??-" as a Trigraph.  See Microsoft VC++ 5.0 help on Trigraphs. BWC.
  //sDate = "??-???-????";
  sDate = "??";
  sDate += "-???";
  sDate += "-????";

  if (NULL != pcTime)
    {
      sDate.SetAt(0, pcTime[ 8]);
      sDate.SetAt(1, pcTime[ 9]);
      sDate.SetAt(3, pcTime[ 4]);
      sDate.SetAt(4, pcTime[ 5]);
      sDate.SetAt(5, pcTime[ 6]);
      sDate.SetAt(7, pcTime[20]);
      sDate.SetAt(8, pcTime[21]);
      sDate.SetAt(9, pcTime[22]);
      sDate.SetAt(10,pcTime[23]);
    }
  return (sDate);
}

int  nDoSystemCommand(const Cstring& rsCommand, int *pnPID)
{
  // Execute an operating system or shell command
  // Wait for return (0 is success, not 0 is failure)
  // REMEMBER: The system routine uses the sh shell.

#if defined(SSI_PC) || defined(WIN32)
  return (system(rsCommand.string()));
#else
  return (system(rsCommand.string()));
/*
  // For Unix systems, try something else
  
  int pid, status;

  if ("" == rsCommand)
    return (-1);
  pid = fork();
  if (-1 == pid)
    return (-1);
  if (0 == pid)
    {
      // It's the child process

      extern char** environ;
      char *argv[4];
      argv[0] = "sh";
      argv[1] = "-c";
      argv[2] = rsCommand.string();
      argv[3] = 0;
      execve("/bin/sh", argv, environ);
      exit(127);
    }

  // Otherwise it's the parent

  if (NULL != pnPID)
    {
      *pnPID = (int)pid;
      cout << "nDo...child pid is " << pid << endl << flush;
    }
  return (0);
*/
#endif
}

#ifdef SSI_PC

int  nDoSystemCommand2(const Cstring& rsCommand, const Cstring& rsOutFile)
{
  // Function does same thing as nDoSystemCommand except that it does not automatically display a cmd.exe 
  // generated DOS console window


        STARTUPINFO si;
        PROCESS_INFORMATION pi;

        memset(&si, 0, sizeof(si));
        si.cb = sizeof(si);
        si.dwFlags |= STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
        si.wShowWindow = SW_HIDE;

        SECURITY_ATTRIBUTES sa;
        sa.nLength = sizeof(SECURITY_ATTRIBUTES);
        sa.lpSecurityDescriptor = NULL;
        sa.bInheritHandle = TRUE;

        DWORD cf = 0;
        cf |= CREATE_NO_WINDOW;

        HANDLE hStdOut = CreateFile( rsOutFile.string(), GENERIC_READ | GENERIC_WRITE, 
                                                                 FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, CREATE_ALWAYS,
                                                                 FILE_ATTRIBUTE_NORMAL, 0 );

        si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        si.hStdOutput = hStdOut;
        si.hStdError = GetStdHandle(STD_ERROR_HANDLE);

        // Run the executable in a new process
        if( CreateProcess(NULL, (LPTSTR)(LPCTSTR) rsCommand.string(), NULL, NULL, true, cf, 
                                          NULL, NULL, &si, &pi) ) {
                WaitForSingleObject( pi.hProcess, INFINITE );
                CloseHandle(pi.hProcess);
                CloseHandle(pi.hThread);
                CloseHandle(si.hStdOutput);
                return 1;
        }
        else {
                return 0;
        }
}
#endif



int  nFileRename(const Cstring& rsFilenameOld, const Cstring& rsFilenameNew)
{
  // Rename rsFilenameOld to rsFilenameNew

  int nStat;

#if defined(SSI_PC) || defined(WIN32)
  // Make sure the filename we are moving to doesn't exist otherwise
  // "rename" will fail.
  if(bFileExists(sTransSymbol(rsFilenameNew)))
     remove(sTransSymbol(rsFilenameNew));
  nStat = rename(sTransSymbol(rsFilenameOld), sTransSymbol(rsFilenameNew));
//  nStat = nDoSystemCommand("rename " + sTransSymbol(rsFilenameOld)
//                         + " " + sTransSymbol(rsFilenameNew));
#else
  // The Unix side could use 'rename' also.
  nStat = nDoSystemCommand("mv " + sTransSymbol(rsFilenameOld)
                           + " " + sTransSymbol(rsFilenameNew));
#endif

  /// RB 06/30/2006 Since we have had problems in the past with nStat,
  // I am adding an different kind of check to see if we have succeded renaming the file
  
  //if (0 != nStat)
  //  {
  //    nStat = 0;
  //    // ALWAYS RETURN SUCCESS!
  //    // Sometimes we have a failure and I am not sure why!
  //  }

  if( bFileExists(rsFilenameOld) )
      return 1;
  /// end RB

  return 0;
}

int  nFileCopy(const Cstring& rsFilenameOld, const Cstring& rsFilenameNew)
{
    int nStat;



#if defined(SSI_PC) || defined(WIN32)
    if( CopyFile(sTransSymbol(rsFilenameOld),
        sTransSymbol(rsFilenameNew), FALSE) ) {
        nStat = 0;
    }
    else {
        nStat = 1;
    }

    //nStat = nDoSystemCommand("copy \"" + sTransSymbol(rsFilenameOld)
    //    + "\" \"" + sTransSymbol(rsFilenameNew) + "\"");
#else
    // Unix
    nStat = nDoSystemCommand("cp " + sTransSymbol(rsFilenameOld)
        + " " + sTransSymbol(rsFilenameNew));
#endif
    return nStat;
}


bool bFileExists(const Cstring& rsFilename)
{
  // Test whether the file with name rsFilename exists
  // by trying to open it

  bool bIsOpen;

#ifdef WIN95

   struct _stat buffer;
   bIsOpen = (-1 != _stat(rsFilename.string(),&buffer));   // file exists

#else

   struct stat buffer;
   bIsOpen = (-1 != stat(rsFilename.string(),&buffer));  // file exists

  // Below only works if file exists and is readable.
  //ifstream oIn( sTransSymbol(rsFilename).string(), ios::in | ios::nocreate);
  //bIsOpen = oIn.rdbuf()->is_open();
  //if (bIsOpen)
  //  oIn.close();

#endif

  return (bIsOpen);
}

int nFileAppendVersion(const Cstring& rsFilename, const bool bCheckEnv)
{
  // Append a version number to the name of an existing file
  // that has no version number on it. As in:
  //+2008-08-20 JWP
  // Two methods of version numbering will be supported. First the original:
  // Old form (used when environment variable DTREK_FILE_VERSIONING is set to 0) 
  // file.ext.1   oldest
  // file.ext.2
  // file.ext     newest  -> file.ext.3

  // New form (used when environment variable DTREK_FILE_VERSIONING is NOT set to 0 or is unset) 
  // file_1.ext   oldest
  // file_2.ext
  // file.ext     newest  -> file_3.ext

  

  // If bCheckEnv is TRUE, first check the environment variable found
  // in Gs_DTREK_OVERWRITE.
  // If env var is set (even to TRUE or FALSE or whatever),
  // then DO NOT append a version number.

  int nStat = 0;

  if (bCheckEnv && ("" != sGetEnv(Gs_DTREK_OVERWRITE)))
    {
      // Environment variable is set, so do not append

      return (0);
    }

  // If the file does not exist, then return success anyways

  if (!bFileExists(rsFilename))
    return (0);

  // Count previous versions of the file

  int     nVersionNum = 1;
  int     nVersionMethod = 1;
  if ("" != sGetEnv("DTREK_FILE_VERSIONING"))
    {
      nStat = atoi(sGetEnv("DTREK_FILE_VERSIONING").string());
      if (0 == nStat)
        nVersionMethod = 0;  // Only switch if exactly matches 0
      nStat = 0;
    }

  Cstring sFilename;
  do
    {
      if (1 == nVersionMethod)
        {
          // New method
          Cstring sFileExtension;
          Cstring sFileBasename;
	  sFileExtension = sFileGetExtension(rsFilename);
	  
	  if ("" != sFileExtension)
	    {
	      sFileBasename = rsFilename;
	      sFileBasename = sFileBasename.before(sFileExtension);
	      sFileBasename = sFileBasename.substr(0,sFileBasename.length()-1);
	      sFilename = sFileBasename
		+ '_' + (Cstring) nVersionNum + '.' + sFileExtension; 
	    }
	  else
	    {
	      sFilename = rsFilename + '_' + (Cstring) nVersionNum;
	    }
        }
      else
        {
          // Original method is the default
          sFilename = rsFilename + '.' + (Cstring) nVersionNum;
        }
      nVersionNum++;
    }
  while (bFileExists(sFilename));

  // Try to rename the file up to 10 times before really an error
  // This is because we randomly get shell error

  nStat              = 1;
  register int i     = 0;
  while ( (10 > i) && (0 != nStat) )
    {
      nStat = nFileRename(rsFilename, sFilename);
      i++;
    }

  if( 0 != nStat )
    {
      cout << "ERROR renaming file " << rsFilename << " to "
           << sFilename << '\n';
    }
  return (nStat);
}


int nFileCreateDirectory(const Cstring& rsDirName)
{
  // Create the directory rsDirName if it does not already exist
  // Return 0 on success, -1 on directory exists and > 0 on other error
  // No test is made to see if the created/existing directory is writable.

#ifdef WIN32

   struct _stat buffer;

   if(-1 != _stat(rsDirName.string(),&buffer)){   // file exists
           if(0 != (_S_IFDIR&buffer.st_mode))
                   return -1;
      else
                        return 2;
   }
   if(0 != _mkdir(rsDirName.string()))
      return 1;
   else
      return 0;

#else

  if (bFileExists(rsDirName))
    {
      // File exists, is it a directory?
      // If it is a directory, then by definition the file rsDirName/. exists, too.

      if (bFileExists(rsDirName + "/."))
        return (-1);
      else
        return (2);
    }

  // Hmmm, we could use the mkdir system routine, rather than a shell call.
  // In that case, be sure to re-do error returns
  return (nDoSystemCommand("mkdir -p " + sTransSymbol(rsDirName)));

#endif
}


int nFileDelete(const Cstring& rsFilename)
{
  // Delete a file with name rsFilename from disk

  // Perhaps we should use the 'unlink()' system routine?

  if (rsFilename.contains('*') || rsFilename.contains('?') )
    {
      // If the input filename contains wildcards, then remove() won't work

#if defined(SSI_PC) || defined(WIN32)
      return (nDoSystemCommand((Cstring)"del \"" + rsFilename + "\""));
#else
      return (nDoSystemCommand((Cstring)"rm -f " + rsFilename));
#endif
    }
  else
    {
      return (remove(rsFilename.string()));
    }
};

Cstring sFileGetDirectory(const Cstring& rsFilename)
{
  // Get the directory name of the input filename

  // If the input filename does not contain the directory specification,
  //   then return the current working directory
  // If the input filename is a directory and does not end with /, :, or ]
  //   then return the parent directory of the input directory.
  // If the input filename is a directory and does end with /, :, or ]
  //   then return the input filename

  Cstring sDir;

  sDir = sTransSymbol(rsFilename);

  // Look for the first /, ] or : starting from last character in the filename

  register int i;
  for (i = sDir.length()-1; i >= 0; i--)
    {
      if (   ('/' == sDir.GetAt(i)) || (']' == sDir.GetAt(i))
          || (':' == sDir.GetAt(i)) || ('\\' == sDir.GetAt(i)) )
        {
          sDir = sDir.before(i+1);
          i = -2;
        }
    }
  if (-1 == i)
    {
      // No evidence of a leading directory was found,
      // so just use current working directory

      sDir = sGetCWD();
    }
  else if ('/' != sDir.GetAt(0) && '\\' != sDir.GetAt(0) && ':' !=
  sDir.GetAt(1))
    {
      // Directory spec is relative, so make absolute

      sDir = sGetCWD() + sDir;
    }

  return (sDir);
}

Cstring sFileGetExtension(const Cstring& rsFilename)
{
  // Get the extension (typically 3 characters after the last . character
  // of the input filename.

  Cstring sFilename;
  sFilename = sTransSymbol(rsFilename);

  // Look for the FIRST . starting from last character in the filename
  // If a . is not found, return empty string

  register int i;
  for (i = sFilename.length()-1; i >= 0; i--)
    {
      if ('.' == sFilename.GetAt(i))
        {
          return (sFilename.after(i));
        }
    }
  return (Cstring(""));

}
Cstring sFileGetBasename(const Cstring& rsFilename)
{
  // Get the basename of the input filename, that is, remove any
  // directory specification

  Cstring sDir;
  sDir = sTransSymbol(rsFilename);

  // Look for the first /, ] or : starting from last character in the filename

  register int i;
  for (i = sDir.length()-1; i >= 0; i--)
    {
      if (   ('/' == sDir.GetAt(i)) ||  (']' == sDir.GetAt(i))
          || (':' == sDir.GetAt(i)) || ('\\' == sDir.GetAt(i)) )
        {
          return (sDir.after(i));
        }
    }
  return (rsFilename);
}

Cstring sFileGetTempName(const char *dir, const char *pfx)
{
  char    *pcTempFile = NULL;
  Cstring sTempFile;

#ifdef __NUTC__
  // NutCracker does not have tempnam(), only tmpnam().
  pcTempFile = tmpnam(NULL);
#elif defined(_WIN32)
  pcTempFile = _tempnam(dir, pfx);
#else
  // Change this because g++ on linux says to use mkstemp()
  sTempFile = sGetEnv("TMPDIR");
  if ("" == sTempFile)
    {
      sTempFile = *dir;
    }
  else
    {
      if (!bFileExists(sTempFile))
        {
          sTempFile = *dir;
        }
    }
  sTempFile = sTempFile + Cstring('/') + Cstring(pfx) + "XXXXXX";
  pcTempFile = new char [sTempFile.length()+1];
  strcpy(pcTempFile, sTempFile.string());
  int nFD = mkstemp(pcTempFile);
  if (-1 == nFD)
    {
      delete [] pcTempFile;
      pcTempFile = NULL;
    }
  else
    {
      close (nFD);
      nFileDelete(Cstring(pcTempFile));
    }
#endif

  if (NULL == pcTempFile)
    sTempFile = "";
  else
    {
      sTempFile = (Cstring) pcTempFile;
      delete [] pcTempFile;
      pcTempFile = NULL;
    }
  return (sTempFile);
}

DTREK_EXPORT Cstring sFileGetFullPath(Cstring& sRelPath)
{
  char* pcRet;
  Cstring sPathBuf;
  sPathBuf = sRelPath;
#if defined(_WIN32)
  const int nMaxPath = 1000;
  char pcBuf[nMaxPath];
  pcRet = _fullpath(pcBuf,sRelPath.string(),nMaxPath);
  if (pcRet)
    sPathBuf = pcRet;
#else
  char *pcBuf;
  pcBuf = new char [MAXPATHLEN];
  pcRet = realpath(sRelPath.string(), pcBuf);
  if (pcRet)
    sPathBuf = pcBuf;
  delete [] pcBuf;
#endif
  return sPathBuf;
};

Cstring sDtrekGetPrefix(void)
{
  return (sGetEnv("DTREK_PREFIX"));
}

static Cstring Gs_sDtrekModuleName = "d*TREK";

Cstring sDtrekGetModuleName(void)
{
  return (Gs_sDtrekModuleName);
}

void  vDtrekSetModuleName(const Cstring& rsString)
{
  if ("" != rsString)
    Gs_sDtrekModuleName = rsString;
}

void vPrintCopyrightInfo(void)
{
     cout << "\n" << sDtrekGetModuleName() << ":  Copyright (c) 2010 Rigaku\n"
          << D_K_DTREKVersion  
 #ifndef WIN32
          << "\nPlease see the file ${DTREK_ROOT}/ACKNOWLEDGEMENTS "
#else
          << "\nPlease see the file ${DTREK_ROOT}\\ACKNOWLEDGEMENTS "
#endif
          << "\nfor further acknowledgements, copyrights and license information.\n\n";
#ifdef SSI_PC
        
// Some CrystalClear-specific stuff, which is using argc and argv

#endif
}


// Get the byte order of the CPU this image is in

int nGetByteOrderCPU(void)
{
  union
    { INT4  l;
      short int s[2];
    } data;

  data.l = 1;
  if ( data.s[0] == 0 && data.s[1] == 1 )
    return eCPU_big_endian;
  else if ( data.s[0] == 1 && data.s[1] == 0 )
    return eCPU_little_endian;
  else
    return eCPU_other_endian;
}

Cstring sGetByteOrderCPU(void)
{
  int nByte_order = nGetByteOrderCPU();
  if (eCPU_big_endian == nByte_order)
    return ((Cstring)D_K_BigEndian);                   // All lower case for Marty's stuff
  else if (eCPU_little_endian == nByte_order)
    return ((Cstring)D_K_LittleEndian);
  else
    return ((Cstring)D_K_OtherEndian);
}

LONG lFileGetSize(const Cstring& rsFilename)
{
  // Return the size of the named file in bytes (characters)
  // On any error getting the size, return a negative result.

  struct stat tStat;
  int nStat;
  nStat = stat(rsFilename.string(), &tStat);
  if (0 != nStat)
    {
      return ((LONG) nStat);
    }
  else
    {
      return ((LONG) tStat.st_size);
    }
}

int nKillProcess(LONG lPID, int nSignal)
{
  // Interface to unix kill routine.
  // Be recursive about it:
  // Look for any process that has lPID as the parent and kill it, too.
  // Then kill original lPID process.

  Cstring sTemp, sTemp2, sTempFile;
  int nStat;
#define KILL_MAX_WORDS 10
  Cstring asWords[KILL_MAX_WORDS];

#ifdef WIN32
  // Could not find an obvious way to kill a process under
  // Windows NT
  // Well, ok.  I suppose that one could do a signal(KILL)
  nStat = -1;
  return nStat;

#else

  LONG lPPID;
  LONG lPIDsub;
  int nWords;

  sTempFile = sFileGetTempName(".", "DTKIL");
#ifdef LINUX
  sTemp     = "ps el > " + sTempFile;
#elif defined(OSX)
  sTemp     = "ps al > " + sTempFile;
#else
  sTemp     = "ps -af > " + sTempFile;
#endif

  nDoSystemCommand(sTemp);
        
  // Now kill the subprocess.  There is a danger that it has started
  // another subprocess.  However if we put this kill before the ps,
  // then the parent PID of any subprocesses are switched to 1 so we
  // do not detect them.

  //cout << "KILLING " << lPID << endl;
  nStat = kill((pid_t)lPID, nSignal);

  int     nPosPID, nPosPPID;
#if defined(LINUX)
          nPosPID   = 2;
          nPosPPID  = 3;
#elif defined(OSX)
          nPosPID   = 1;
          nPosPPID  = 2;
#else
          nPosPID   = 1;
          nPosPPID  = 2;
#endif
  ifstream oIn(sTempFile);
  if ((oIn.rdbuf()->is_open()))
    {
      while (!oIn.eof())
        {
          getline(oIn, sTemp);
          if (!oIn.eof())
            {
              nWords = split(sTemp, asWords, KILL_MAX_WORDS, " \t\n");
              if ((nPosPPID+1) <= nWords)
                {

                  lPPID   = (LONG)atoi(asWords[nPosPPID].string());
                  lPIDsub = (LONG)atoi(asWords[nPosPID].string());
                  //cout << "PIDsub, parent " << lPIDsub << ' ' << lPPID << endl;
                  if (lPID == lPPID)
                    {
                      //cout << "SubPID is: " << lPIDsub << endl;

                      // Notice recursion!!!

                      nKillProcess(lPIDsub, nSignal);
                    }
                }
            }
        }
      oIn.close();
      nFileDelete(sTempFile);
    }

  return (nStat);

#endif   /* #ifdef WIN32 */
}

int  nDoSystemCommand(const Cstring& rsCommand, FILE **ppSubprocessStdIn, int *pnPID)
{
  // Run a command in a subshell, return a file handle to the standard input
  // of the subshell
  // Note that *pSubprocessStdIn MUST be closed with pclose().

#ifdef WIN32
  *ppSubprocessStdIn = _popen(rsCommand.string(), "w");
#else
  *ppSubprocessStdIn = popen(rsCommand.string(), "w");
#endif

  if (NULL == *ppSubprocessStdIn)
    {
      cout << "ERROR in nDoSystemCommand popen failed!\n";
      return (-1);
    }
  else
    {
      return (0);
    }

}

/****************************************************************************
 * Returns the time of day in seconds since midnight UTC 01 Jan 1970.       *
 ****************************************************************************/
double dGetTimeOfDay(void)
{
   double dTime;

#ifdef WIN32

   struct _timeb oTime;
   _ftime( &oTime );
   dTime = oTime.time+oTime.millitm/1000.0;

#else
#ifndef LINUX
   struct timeval oTime;
#if defined(OSF1) || defined(SUNOS) || defined(OSX)
   struct timezone TZ;
   gettimeofday( &oTime, &TZ);
#else
   gettimeofday( &oTime );
#endif
   dTime = oTime.tv_sec+oTime.tv_usec/1000000.0;
#else

   struct timeb oTime;
   ftime( &oTime );
   dTime = oTime.time+oTime.millitm/1000.0;

#endif

#endif

   return dTime;
}


Cstring sFileBuildName(Cstring &rsDirectory,
                             Cstring &rsFilename)
{
   int nLen;
   Cstring sResult;

   sResult = sTransSymbol(rsDirectory);
   nLen = sResult.length();

#ifdef WIN32
   if('\\' != sResult.GetAt(nLen-1))
      sResult += '\\';
#else
   if('/' != sResult.GetAt(nLen-1))
      sResult += '/';
#endif

   sResult += rsFilename;

   return sResult;
}

Cstring sFileBuildName(char *pcDirectory,
                             Cstring &rsFilename)
{
   Cstring s = pcDirectory;
   return sFileBuildName(s,rsFilename);
}

Cstring sFileBuildName(Cstring &rsDirectory,
                             char *pcFilename)
{
   Cstring s = pcFilename;
   return sFileBuildName(rsDirectory,s);
}

Cstring sFileBuildName(char *pcDirectory,
                             char *pcFilename)
{
   Cstring s = pcDirectory;
   Cstring ss = pcFilename;
   return sFileBuildName(s,ss);
}

double dGetFreeSpaceOnDisk(Cstring &sFileOnDisk)
{
   int nStatus;
   int n;
   double dFreeBytes;
   Cstring sFile;

   if("" == sFileOnDisk){
      sFile = sGetCWD();
   }
   else
      sFile = sFileOnDisk;

#ifdef WIN32

   char c;
   Cstring sDrive;

   // Make sure that the file does not have a directory slash at the end of
   // the name.  stat() will fail if it does.   Also if just a drive letter,
   // with a colon, then tack on a backslash to keep stat happy.

   n = sFile.length()-1;
   if('\\' == sFile.GetAt(n))
       sFile += '.';
   else if(':' == sFile.GetAt(n))
       sFile += '\\';

// Find out which drive we are on.

   struct stat oStat;
   nStatus = stat(sFile.string(),&oStat);
   if(0 != nStatus)
       return 0.0;
   c = oStat.st_dev+'A';
   sDrive = c;
   sDrive += ":\\";

   // Get free disk space info for that drive

   typedef BOOL (__stdcall *EX) (LPCTSTR,PULARGE_INTEGER,PULARGE_INTEGER,PULARGE_INTEGER);

   EX pGetDiskFreeSpaceEx = (EX) GetProcAddress(GetModuleHandle("kernel32.dll"),
       "GetDiskFreeSpaceExA");

   if(NULL != pGetDiskFreeSpaceEx){
       __int64 TotalBytesAvailableToUser,TotalNumberOfBytes,TotalNumberOfFreeBytes;
       BOOL b = pGetDiskFreeSpaceEx(sDrive.string(),
           (PULARGE_INTEGER)&TotalBytesAvailableToUser,
           (PULARGE_INTEGER)&TotalNumberOfBytes,
           (PULARGE_INTEGER)&TotalNumberOfFreeBytes);
       if(b)
           dFreeBytes = (double)TotalBytesAvailableToUser;
       else
           dFreeBytes = 0.0;
   }
   else{
       DWORD dwSectorsPerCluster,dwBytesPerSector,dwFreeClusters,dwTotalClusters;
       if(FALSE == GetDiskFreeSpace(sDrive.string(),&dwSectorsPerCluster,
           &dwBytesPerSector,&dwFreeClusters,&dwTotalClusters))
           dFreeBytes = 0.0;
       else{
           dFreeBytes = (double)dwFreeClusters;
           dFreeBytes *= (double)dwSectorsPerCluster;
           dFreeBytes *= (double)dwBytesPerSector;
       }
   }

   /*char c;
   struct stat oStat;
   DWORD dwSectorsPerCluster,dwBytesPerSector,dwFreeClusters,dwTotalClusters;
   Cstring sDrive;

// Make sure that the file does not have a directory slash at the end of
// the name.  stat() will fail if it does.   Also if just a drive letter,
// with a colon, then tack on a backslash to keep stat happy.

   n = sFile.length()-1;
   if('\\' == sFile.GetAt(n))
      sFile += '.';
   else if(':' == sFile.GetAt(n))
      sFile += '\\';

// Find out which drive we are on.

   nStatus = stat(sFile.string(),&oStat);
   if(0 != nStatus)
      return 0.0;
   c = oStat.st_dev+'A';
   sDrive = c;
   sDrive += ":\\";

// Get free disk space info for that drive

   if(FALSE == GetDiskFreeSpace(sDrive.string(),&dwSectorsPerCluster,
   &dwBytesPerSector,&dwFreeClusters,&dwTotalClusters))
      return 0.0;

// Piece together the number of free bytes.

   dFreeBytes = (double)dwFreeClusters;
   dFreeBytes *= (double)dwSectorsPerCluster;
   dFreeBytes *= (double)dwBytesPerSector;*/

#else

   struct statfs oStatfs;

#if defined(LINUX) || defined(OSX)
   nStatus = statfs(sFile.string(),&oStatfs);
#else
   nStatus = statfs(sFile.string(),&oStatfs,sizeof(struct statfs),0);
#endif

   if(0 != nStatus)
      return 0.0;

// statfs returns info in blocks, so convert to bytes.


#ifdef OSF1
   dFreeBytes = (double)oStatfs.f_bfree*(double)oStatfs.f_fsize;
#elif defined(SUNOS)
   dFreeBytes = (double)oStatfs.f_bfree*512.0;
#else
   dFreeBytes = (double)oStatfs.f_bfree*(double)oStatfs.f_bsize;
#endif

#endif   /* ifdef(WIN32) */

   return dFreeBytes;
}
#if defined(__cdecl)
double __cdecl dGetFreeSpaceOnDisk(const char *pcFileOnDisk)
#else
double dGetFreeSpaceOnDisk(const char *pcFileOnDisk)
#endif
{
   Cstring s = pcFileOnDisk;
   return dGetFreeSpaceOnDisk(s);
}

#if defined(WIN32)

int nint(double d)
{
   int i;
   if(0 > d)
      i = -((int)(-d+0.5));
   else if(0 < d)
      i = (int)(d+0.5);
   else
      i = 0;
   return i;
}
/*
long rint(double d)
{
   long i;
   if(0 > d)
      i = -((long)(-d+0.5));
   else if(0 < d)
      i = (long)(d+0.5);
   else
      i = 0;
   return i;
}
*/
#endif

/*
Cstring GetInstallDirectory(void)
{
        Cstring csRet = ".";

#ifdef SSI_PC
        // Look in the registry for the location of the configuration file.
        TCHAR SUBKEYNAME[] = "SOFTWARE\\Rigaku MSC\\CrystalClear\\1.00.000\\Location of Files";
        HKEY hkResult;
        DWORD vType;
        DWORD dSize = _MAX_PATH;
        BYTE Path[_MAX_PATH];
        Path[0]= '\0';

        long ret = RegOpenKeyEx(HKEY_LOCAL_MACHINE,     SUBKEYNAME,
                                                        0, KEY_QUERY_VALUE,     &hkResult);

        if(ret == ERROR_SUCCESS)
        {
                if(hkResult)
                {
                                // Get the path of the install directory.
                                RegQueryValueEx(hkResult, "Install Directory",
                                                                NULL, &vType, (LPBYTE)Path, &dSize);
                }
                RegCloseKey(hkResult);
        }
        csRet = (char*)Path;
#endif

        return csRet;
}
*/

Cstring
sBuildScanTemplate(const Cstring &rsFilename, const int nPlaces, 
                   int *pnSeqNum)
{
  // Build and return a d*TREK scan template from the input filename
  // If pnSeqNum is not NULL, then return the sequence number
  // of rsFilename, too.

  int     nTemp, i, j;
  Cstring sScanTemplate;

  if ("" == rsFilename)
    return ("");

  sScanTemplate = rsFilename;

  // Change last nPlaces digits in rsFilename to ?s,
  // but treat MAR IP images special
  // since they usually end in .mar???? or .pck???? where ???? is the image
  // size in pixels.
  // Also treat images ending in .bz2 as special.

  int nLast;
  int nSeqNum = 0;
  int nTen    = 1;

  nLast = sScanTemplate.find(".mar") - 1;
  if (1 > nLast)
    nLast = sScanTemplate.find(".MAR") - 1;
  if (1 > nLast)
    nLast = sScanTemplate.find(".pck") - 1;
  if (1 > nLast)
    nLast = sScanTemplate.find(".PCK") - 1;
  if (1 > nLast)
    {
      nLast = sScanTemplate.find(".bz2") - 1;
      //cout << "nLast w/bz2 is: " << nLast << endl;
      //cout << "scantemp is: " << sScanTemplate << endl;
    }

  if (sScanTemplate.contains(".bz2"))
    {
      if (nLast != (sScanTemplate.length() - 5))
        {
          nLast = sScanTemplate.length()-1;
        }
    }
  else if ( (1 > nLast) || (nLast != (sScanTemplate.length() - 9)) )
    {
      nLast = sScanTemplate.length()-1;
    }

  nTemp = 0;
  bool bFirstDashFound = FALSE;

  for (i = nLast; (nTemp < nPlaces) && (i >= 0); i--)
    {
      j = sScanTemplate.GetAt(i);
      //cout << "i is: " << i << " and j is " << j << endl;
      if ( (j >= '0') && (j <= '9') )
        {
          nSeqNum = nSeqNum + ((int) (j - char('0')) * nTen);
          nTen   *= 10;
          sScanTemplate.SetAt(i, '?');
          nTemp++;
        }
      else if ( ('a' <= j) && ('z' >= j) && (0 < nTemp) && (3 >= nPlaces) )
        {
          // This should probably only be used for BS images
          // and not for templates with nPlaces > 3.
          // A letter between a and z, inclusive,
          // but only after an earlier digit was found

          nSeqNum   = nSeqNum + nTen * (j - (int) 'a' + 10);
          nTen     *= 10;
          nTemp++;
        }
      else if (  ('-' == j) || ('_' == j) )
        {
//+2011-08-17 JWP
          // Stop looking for digits if we get to a 2nd underscore or dash
	  // and we have not finished nPlaces yet
          // (we may wish to change this in the future)

          if (bFirstDashFound) i = -1;
	  bFirstDashFound = TRUE;
//-2011-08-17 JWP
        }  
      else if (  ('/' == j) || ('\\' == j) || (']' == j) || (':' == j) )
        {
          // Stop looking for digits if we get to a directory indicator
          // (we may wish to change this in the future)

          i = -1;
        }
    }

  // Prepend a directory specification if it is not already there

  if (sScanTemplate == sFileGetBasename(sScanTemplate))
    sScanTemplate = sGetCWD() + sScanTemplate;

  if (NULL != pnSeqNum)
    *pnSeqNum = nSeqNum;
  return (sScanTemplate);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring sBuildScanPrefix(const Cstring& sPre,int nScanNumber)
{
    Cstring     sEmpty("");
    
    return sBuildStrategyScanPrefix(sEmpty, sPre, nScanNumber);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring sBuildStrategyScanPrefix(const Cstring& sStrategyPrefix, const Cstring& sPre, int nScanNumber)
{
    Cstring   sPrefix(sStrategyPrefix);

    if( "" == sStrategyPrefix && 0 < nScanNumber || "" != sStrategyPrefix )
    { 
        char      pcTemp[20];
        sprintf(pcTemp,"S%d_", nScanNumber);
        sPrefix += pcTemp;
    }

    sPrefix += sPre; 

    return sPrefix;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

Cstring g_sModifiedCommandLine;
itr<char*> g_apcModifiedCommandLine;


int nModifyCommandLine(unsigned int& nArgc,char*** pargv)
{
    Cstring         sTemp;
    Cstring         sArg;
    Cstring         sFile;
    int             nLine;
    bool            bRecording;
    std::vector<Cstring>    asModifiedCommandLine;
    char**          argv = *pargv;
    const int       nMaxBufLength = 200;
    char            pcBuf[nMaxBufLength + 1];
    int             nPrevLength;
    FILE*           pFLog;
    char*           pcPoint;
    int             nArg;
    bool            bHadContinuationMarker;

    if ((nArgc==2) && ((sFile = argv[1]).contains(".log"))) {
        // User has requested to run the program from a log file
        pFLog = fopen(sFile,"rt");
        if (!pFLog)
            return 0;
        // We don't set the recording flag until we find the "comand line:" string.
        bRecording = false;
        nLine = 0;        
        while ((!feof(pFLog)) && (nLine<200)) {
            bHadContinuationMarker = false;
            if (fgets(pcBuf,nMaxBufLength,pFLog)) {
                sTemp = pcBuf;
                if (bRecording) {
                    nPrevLength = asModifiedCommandLine.size();
                    while (sTemp.length()) {
                        sTemp = sTemp.after(' ');
                        sArg = sTemp.before(' ');
                        sArg = sArg.before('\n');
                        if (sArg == "\\")
                            bHadContinuationMarker = true;
                        if ((!sArg.length()) || (sArg == "\\"))
                            continue;
                        else
                            asModifiedCommandLine.push_back(sArg);
                    };
                    // If we didn't add anything, then we are done.
                    if ((asModifiedCommandLine.size()==nPrevLength) || (!bHadContinuationMarker)) {
                        break;
                    };
                } else if (sTemp.contains("Command line:"))
                    bRecording = true;

                nLine++;               
            } else
                continue;
        };
        fclose(pFLog);
        if (asModifiedCommandLine.size()) {
            

            printf("INFO:  Log file '%s' is for d*TREK module '%s'!\n",sFile.string(),asModifiedCommandLine[0].string());
            g_sModifiedCommandLine = "";
            -g_apcModifiedCommandLine;
            nArgc = asModifiedCommandLine.size();
            for (nArg = 0; nArg < nArgc; nArg++) {
                g_sModifiedCommandLine += asModifiedCommandLine[nArg];
                g_sModifiedCommandLine += '\x0';                    
            };
            for (nArg = 0,pcPoint = g_sModifiedCommandLine.string();nArg <nArgc; nArg++) {
                g_apcModifiedCommandLine + (pcPoint);
                pcPoint += (strlen(pcPoint) + 1);
            };
            *pargv = &g_apcModifiedCommandLine[0];

            // Write out the modified command line.
            FILE* pFOut;
            pFOut = fopen("logcmd","wt");
            for (nArg = 0; nArg < nArgc; nArg++) {
                fprintf(pFOut,"%s ",asModifiedCommandLine[nArg].string());
            };
            fprintf(pFOut,"\n");
            fclose(pFOut);
            

        };


    };
    return 0;
};

Cstring g_sCommandLine;

Cstring sGetCommandLine(int argc, char *argv[], const int nMaxWidth)
{
    // Get command line arguments as a Cstring for possible output
    int     i;
    int     nLength;
    Cstring sTemp,sArg;

    g_sCommandLine = "";

    // Copy over the command line (without the '\' characters)
    for(i = 0; i < argc;i++)
    {
        g_sCommandLine += argv[i];
        
        g_sCommandLine += " ";
    }

    // Copy command line files to sRefineFiles
    if( 0 < argc ) 
        sArg = (Cstring)argv[0];

    nLength  = sArg.length();

    for (i = 1; (i < argc); i++)
    {
        sTemp = Cstring(' ') + (const char*) argv[i];
        
        if( 0 != nMaxWidth && (sTemp.length() + nLength) > nMaxWidth )
        {
            // String line exceeded maximum allowed width so insert newline
            sArg  += " \\\n   ";
            nLength = 3;
        }
        
        sArg = sArg + sTemp;
        
        nLength += sTemp.length();
    }
    
    //  cout << sArg << endl;
    return sArg;
}

Cstring sTabFormat(Cstring& sString, int nTabCount,int nMaxWidth)
{
    static Cstring sBuf;
    static Cstring sTemp;
    int nTab;
    int nWidth;
    int nLastws;
    int nLastwsOut;
    char cc;
    int nx;

    // Macro to add the tabs when we need them.
#define ADDTABS for (nTab=0;nTab<nTabCount;nTab++) { sBuf += " "; nWidth++; };


    sBuf = "";
    nWidth =0;
    nLastws = -1;
    
    ADDTABS;
    for (nx=0;nx<sString.length();nx++) {
      //cc = sString.string()[nx];
      cc = sString.GetAt(nx);
        if (cc=='\n') {
            sBuf +="\n";
            nWidth = 0;
            nLastws = -1;
            ADDTABS;
        } else if (cc==' ') {
            nLastws = nx;
            nLastwsOut = sBuf.length();
            sBuf += " ";
            nWidth ++;
        } else if (nWidth>nMaxWidth) {
            if (nLastws!= -1) {

                sTemp="";
                sTemp.append(sBuf,0,nLastwsOut);
                sBuf = sTemp;
                sBuf+="\n";
                nx = nLastws;
                nWidth=0;
                nLastws = -1;
                ADDTABS;
            } else {
                sBuf += cc;
                nWidth ++;
            };
        } else {
            sBuf += cc;
            nWidth ++;
        };
    };
                    
#undef ADDTABS
    return sBuf;
};


float
fSwapFloat(float fFloat)
{
  union _tagSwap
  {
    float fTemp;
    UINT4 ulTemp;
  } tSwap;
  
  tSwap.fTemp = fFloat;
  tSwap.ulTemp = (tSwap.ulTemp << 24) | ((tSwap.ulTemp << 8) & 0x00ff0000) |
                  ((tSwap.ulTemp >> 8) & 0x0000ff00) | (tSwap.ulTemp >> 24);

  return tSwap.fTemp;
}

float
fVAXtoIEEE(float *pfTemp)           // Swaps and converts?
{
  union {
    float    fFloat;
    char     cByte[4];
  } tData;
  char cTemp;

  // This algorithm comes from some R-AXIS FORTRAN software.
  if (*pfTemp == 0.0) return (0.0);

  tData.fFloat   = *pfTemp;
  cTemp          = tData.cByte[0];
  tData.cByte[0] = tData.cByte[1] - 1;
  tData.cByte[1] = cTemp;
  cTemp          = tData.cByte[3];
  tData.cByte[3] = tData.cByte[2];
  tData.cByte[2] = cTemp;
#if defined(WIN32) || defined(LINUX) || defined(OSF1)
  // Must swap float after converting VAX to IEEE on a IntelPC
  tData.fFloat   = fSwapFloat(tData.fFloat);
#endif

  return ( tData.fFloat );
}

float
fIEEEtoVAX(float *pfTemp)           // Swaps and converts?
{
  union {
    float    fFloat;
    char cByte[4];
  } tData;
  char cTemp;

  if (*pfTemp == 0.0) return (0.0);

  // This algorithm comes from some R-AXIS FORTRAN software.
  tData.fFloat   = *pfTemp;
  cTemp          = tData.cByte[0];
  tData.cByte[0] = tData.cByte[1] + 1;
  tData.cByte[1] = cTemp;
  cTemp          = tData.cByte[3];
  tData.cByte[3] = tData.cByte[2];
  tData.cByte[2] = cTemp;

  return ( tData.fFloat );
}


int nGetBatchNumber(Cstring& sBatch, Cstring& sScanID)
{
    sScanID = "";
    
    int         nScanID = 0;

    int         nBatchNameLength = sBatch.length();
    int         nNumDigits = 0;

    Cstring     sTemp(sBatch);
    for(int nx=nBatchNameLength-1; nx >= 0; nx--)
    {
        sTemp.upcase();

        if( sTemp.GetAt(nx) > '9' ||  sTemp.GetAt(nx) < '0' ) 
            break;
        
        nNumDigits++;
    }

    if( nNumDigits == 0 )
        nScanID = 0;
    else
    {
        // Parse up to last three digits.
        nNumDigits = min(4, nNumDigits);
        if( 1 != sscanf(sBatch.substr(nBatchNameLength - nNumDigits, nNumDigits).string(),"%d",&nScanID)) 
            return -1;
    }
    
    sScanID = sBatch.substr(0, nBatchNameLength - nNumDigits);
    
    return  nScanID;    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void vMakeBatchName(int nBatchNumber, const char* pcScanPrefix, Cstring& sBatchName)
{
    char        szBuffer[256];
    sprintf(szBuffer, "%s%04d", pcScanPrefix, nBatchNumber);

    sBatchName = szBuffer;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


Cstring* g_psPrintfLog = NULL;
int lprintf(const char* pcFormat,...) {
    char pcBuf[1000];  // Hope we won't need a buffer bigger than this!
    int nStat;
    va_list oArgList;
    
    va_start(oArgList,pcFormat);    
    
    nStat = vsprintf(pcBuf,pcFormat,oArgList);
    printf("%s",pcBuf);
    if (g_psPrintfLog)
        (*g_psPrintfLog) += pcBuf;
    return nStat;
};


DTREK_EXPORT int nWait(int nMillSecond) {
#ifdef WIN32
    Sleep(nMillSecond);
#else
    Cstring sTemp;
    sTemp = "sleep ";
    sTemp += (nMillSecond/1000);
    nDoSystemCommand(sTemp);
#endif
    return 0;
};


int nPrintTable(Cstring& sTableTitle,
                std::vector<Cstring>& asTableHeaders,
                std::vector<Cstring>& asTableFormat,
                itr<double>& afValues,
                bool bLastRowIsSummary)
{
    int nTotalTableLength;
    int nCol;
    int nDataCols;
    int nx,ny,nz;
    Cstring sTemp;
    double f0;
    
    itr<int> anHeaderCol;
    itr<int> anType;
    std::vector<Cstring> asFragment;
    const int nChar = 1;
    const int nFloat = 2;
    const int nInt = 3;
    const int nNothing = 4;
    
    nTotalTableLength = 0;
    nDataCols = 0;
    for (nCol = 0; nCol < asTableHeaders.size();nCol++) {
        int nType;
        int nLastFragmentEnd;
        int nThisFragmentEnd;
        
        nTotalTableLength += asTableHeaders[nCol].before('\n').length();
        
        // Must discover how many numerical collumns we have.  This is very important,
        // since a single heading might have more than one data item.
        
        nLastFragmentEnd = -1;
        sTemp = asTableFormat[nCol];
        while (sTemp.contains('%')) {
            nThisFragmentEnd = sTemp.find('%',sTemp.find('%')+1);
            nType = nNothing;
            ny = sTemp.find('c');
            if ((nThisFragmentEnd==-1) || ((ny < nThisFragmentEnd) && (ny != -1))) {
                nType = nChar;
                nThisFragmentEnd = ny;
            };
            ny = sTemp.find('f');
            if ((nThisFragmentEnd==-1) || ((ny < nThisFragmentEnd) && (ny != -1))) {
                nType = nFloat;
                nThisFragmentEnd = ny;
            };
            ny = sTemp.find('d');
            if ((nThisFragmentEnd==-1) || ((ny < nThisFragmentEnd) && (ny != -1))) {
                nType = nInt;
                nThisFragmentEnd = ny;
            };
            if (nThisFragmentEnd == -1)
                return 1;
            else {
                sTemp = sTemp.after(nThisFragmentEnd);
                nThisFragmentEnd += (nLastFragmentEnd + 1);
            };

            anHeaderCol + nCol;
            anType + nType;                
            
            asFragment.push_back(asTableFormat[nCol].before(nThisFragmentEnd + 1).after(nLastFragmentEnd));
            nLastFragmentEnd = nThisFragmentEnd;
            if (nType != nNothing)
                nDataCols ++;
        };
        anType + nNothing;
        anHeaderCol + nCol;
        asFragment.push_back(asTableFormat[nCol].after(nLastFragmentEnd));
        
    };
    // The # of columns must be a multiple of afValues.
    if ((!anHeaderCol.size()) || (afValues.size() % nDataCols))
        return 1;

    printf("%s\n",sTableTitle.string());
    for (nx = 0; nx < nTotalTableLength; nx++)
        printf("-");
    printf("\n");
    for (nx = 0; nx < 2; nx++) {
        for (nCol = 0; nCol < asTableHeaders.size();nCol++) {
            int nFirstRowLength = asTableHeaders[nCol].before('\n').length();
            int nSecondRowLength= asTableHeaders[nCol].after('\n').length();

            if (nx == 0)
                printf("%s",asTableHeaders[nCol].before('\n').string());
            else {                    
                printf("%s",asTableHeaders[nCol].after('\n').string());
                while (nSecondRowLength < nFirstRowLength) {
                    printf(" ");
                    nSecondRowLength++;
                };
            };
        };
        printf("\n");
    };
    for (nx = 0; nx < nTotalTableLength; nx++)
        printf("-");
    printf("\n");

    // Now print the table.
    
    int nDataPoint;

    nDataPoint = 0;
    while (nDataPoint < afValues.size()) {
        for (nCol = 0; nCol < asTableHeaders.size();nCol++) {
            for (nx = 0; nx < anHeaderCol.size(); nx++) {
                const double g_fTableAverage = 1.123456e-90;

                if (nDataPoint < afValues.size()) {
                    if (afValues[nDataPoint] == g_fTableAverage) {
                        f0 = 0.0;
                        for (ny = nDataPoint - nDataCols, nz = 0; ny >= 0; ny -= nDataCols,nz++) {
                            f0 += afValues[ny];
                        };
                        f0/=max(1,nz);
                    } else if (afValues[nDataPoint] == g_fTableMax) {
                        for (ny = nDataPoint - nDataCols, nz = 0; ny >= 0; ny -= nDataCols,nz++) {
                            if (nz == 0)
                                f0 = afValues[ny];
                            else
                                f0 = max(afValues[ny],f0);
                        };
                    } else if (afValues[nDataPoint] == g_fTableMin) {
                        for (ny = nDataPoint - nDataCols, nz = 0; ny >= 0; ny -= nDataCols,nz++) {
                            if (nz == 0)
                                f0 = afValues[ny];
                            else
                                f0 = min(afValues[ny],f0);
                        };
                    } else if (afValues[nDataPoint] == g_fTableSum) {
                        f0 = 0.0;
                        for (ny = nDataPoint - nDataCols, nz = 0; ny >= 0; ny -= nDataCols,nz++) {
                            f0 += afValues[ny];
                        };
                    } else
                        f0 = afValues[nDataPoint];
                };


                if (anHeaderCol[nx] == nCol) {
                    if (anType[nx] == nChar) {
                        printf(asFragment[nx].string(),(char) f0);
                        nDataPoint++;
                    } else if (anType[nx] == nFloat) {
                        printf(asFragment[nx].string(),(double) f0);
                        nDataPoint++;
                    } else if (anType[nx] == nInt) {
                        printf(asFragment[nx].string(),(int) f0);
                        nDataPoint++;
                    } else if (anType[nx] == nNothing) {
                        printf(asFragment[nx].string());
                    };
                };
            };
        };
        printf("\n");
        if (nDataPoint + nDataCols == afValues.size()) {
            for (nx = 0; nx < nTotalTableLength; nx++)
                printf("-");
            printf("\n");
        };
    };
    for (nx = 0; nx < nTotalTableLength; nx++)
        printf("-");
    printf("\n");

    return 0;
};




int nTokenize(Cstring& sStringIn, std::vector<Cstring>& asTokensOut, std::vector<Cstring>* pasSymbols)
{
    char* pcToks = "\t\n\r ";
    char* pcTok;
    Cstring sTemp;

    asTokensOut.clear();
    for (pcTok = strtok(sStringIn.string(),pcToks);pcTok;pcTok = strtok(NULL,pcToks)) {
        // Further refine the token.
        if (pasSymbols) {
            int nMinIndex;
            int nMinSymbol;
            int nSymbol;
            int nIndex;
            Cstring sRemainingToken;

            sRemainingToken = ((Cstring) pcTok);
            
            while (sRemainingToken.length()) {
                nMinIndex = -1;
                for (nSymbol = 0; nSymbol < pasSymbols->size(); nSymbol++) {
                    if ((nIndex = sRemainingToken.find((*pasSymbols)[nSymbol]))>=0) {
                        if ((nMinIndex == -1) || (nIndex < nMinIndex)) {
                            nMinIndex = nIndex;
                            nMinSymbol = nSymbol;
                        };
                    };
                };
                if (nMinIndex == -1) {
                    asTokensOut.push_back(sRemainingToken);
                    sRemainingToken = "";
                } else {
                    sTemp = sRemainingToken.before(nMinIndex);
                    if (sTemp.length())
                        asTokensOut.push_back(sTemp);
                    asTokensOut.push_back((*pasSymbols)[nMinSymbol]);
                    sRemainingToken = sRemainingToken.after((*pasSymbols)[nMinSymbol]);
                };
            };
        } else
            asTokensOut.push_back((Cstring)pcTok);
    };
    return 0;
};

//
void vModifyDISPLAY(void)
{
  // This routine is supposed to examine the X Windows related DISPLAY
  // environment variable and change it to be the IP number and not the
  // alphabet hostname.  For example, mydisplay.mynode.com:0.0 is changed to
  // 192.192.1.5:0.0
  // If the DISPLAY env variable is already all-numeric, 
  // or if the DISPLAY env variable is empty,
  // or if the DISPLAY env variable has no characters before the colon,
  // or if the DISPLAY env variable has no colon,
  // then nothing is done.
  // This routine creates a temporary file and deletes it.

  //+JWP 2006-08-04
  // New addition, check if XKEYSYMDB environment variable is set and valid.

  Cstring sXkeysymDB;

  sXkeysymDB = sGetEnv("XKEYSYMDB");
  if ("" == sXkeysymDB)
    {
      // Not set, try to set it if XKeysymDB not found in /usr/X11R6/lib/X11

      sXkeysymDB = "/usr/X11R6/lib/X11/XKeysymDB";
      if (!bFileExists(sXkeysymDB))
        {
          // File does not exist in the usual place, so look in /usr/share...

          sXkeysymDB = "/usr/share/X11/XKeysymDB";
          if (bFileExists(sXkeysymDB))
            nPutEnv("XKEYSYMDB", sXkeysymDB);
        }
    }
  //-JWP 2006-08-04
#ifndef OSX
  // Do not try this if Mac/OSX version

  Cstring sDisplay, sBeforeColon, sAfterColon, sTempFile, sCommand;
  sDisplay     = sGetEnv("DISPLAY");

  sAfterColon  = sDisplay.after(':');
  sBeforeColon = sDisplay.before(':');

  if (0 >= sDisplay.find(':'))
    return;

  // See if sBeforeColon has any characters besides . 0-9.

  int nx;
  bool bIPnumOnly = TRUE;
  for (nx = 0; nx < sBeforeColon.length(); nx++)
    {
      if (    ('.' > sBeforeColon.GetAt(nx))
           || ('9' < sBeforeColon.GetAt(nx))) 
        bIPnumOnly = FALSE;
    }
  if (bIPnumOnly)
    {
      return;
    }
  else
    {
      sTempFile = sFileGetTempName(".", "TMPDPY");
      // What if "host" command does not exist?
      Cstring sHostCommand = "DT_HOSTCOMMAND";
      sCommand = sGetEnv(sHostCommand); 
      if ("" == sCommand)
	sCommand = "host";
      else if ("RETURN" == sCommand)
	return;
      sCommand = sCommand + " " + sBeforeColon + " 1> " + sTempFile + " 2>&1 ";
      //cout << "HOSTCOMMAND: >>" << sCommand << "<<\n";
      nDoSystemCommand(sCommand);
      FILE    *fp;
      char    cLine[85];
      fp = fopen(sTempFile.string(), "r");
      if (fp)   
        {
          fgets(cLine, 80, fp);
          //cout << "HOSTRESPONSE: >>" << cLine << "<<\n";
          sBeforeColon = Cstring(cLine);

          // What if output of "host" command does not have "address" in it?
          nx = 0;
          if (0 >= sBeforeColon.length())
            {
              // String was empty
              sBeforeColon = "localhost has address 127.0.0.1\n";
              printf("WARNING setting display IP to 127.0.0.1.\n");
            }
          else if (sBeforeColon.contains("not found"))
            {
              if (sBeforeColon.contains("local"))
                {
                  // Host command returned error, try setting to 127.0.0.1
                  printf("WARNING setting DISPLAY to 127.0.0.1:%s\n",
                         sAfterColon.string());
                  sBeforeColon = "address 127.0.0.1\n";
                }
              else
                {
                  printf("ERROR setting DISPLAY, left unchanged as %s\n",
                         sDisplay.string());
                  nx = 1;
                }
            }
          if (0 == nx)
            {
              sBeforeColon = sBeforeColon.after("address ");
              sBeforeColon = sBeforeColon.before('\n');
              if (6 <= sBeforeColon.length())
                {
                  // Change DISPLAY environment variable
                  sDisplay = sBeforeColon + ':' + sAfterColon;
                  //cout << "sDisplay>>" << sDisplay << "<<\n";

                  //***********************************************************
                  // This is where the environment variable is actually changed
                  //***********************************************************

                  nPutEnv("DISPLAY", sDisplay);
                }
            }
        }
      fclose (fp);
      (void) nFileDelete(sTempFile);
    }
#endif
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef VC9
char* dt_strcpy(char *strDestination, const char *strSource)
{
        if( 0 != strcpy_s(strDestination, sizeof(strSource), strSource) )
                return NULL;

        return strDestination;
}
///////////////////////////////////////////////////////////////////////////////////////////
char* dt_strncpy(char *strDest, const char *strSource, size_t count)
{
        if( 0 != strncpy_s(strDest, sizeof(strDest), strSource, count) )
                return NULL;

        return strDest;
}
///////////////////////////////////////////////////////////////////////////////////////////
char* dt_strtok(char *strToken, const char *strDelimit)
{
        static char*    strContext = NULL;
                
        return strtok_s(strToken, strDelimit, &strContext);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
// Below are some function for managing debug output
// 1. Specify the debug file path - call set_debug_file_path()
// 2. Specify (optionally) an environment variable name to be tested before each output - call set_debug_env_var_name()
// 3. Create debug file - call create_debug_file()
// 4. Print debug information into the file - call print_debug_string()
// 5. Close file - call close_debug_file()
// 6. Delete file - call delete_debug_file()
static  FILE*   s_fpDebugFile = NULL;
static  Cstring s_sDebugFilePath("");
static  Cstring s_sDebugEnvVarName("");

void set_debug_file_path(Cstring& sPath)
{
   s_sDebugFilePath = sPath;
}

void set_debug_env_var_name(Cstring& sVar)
{
    s_sDebugEnvVarName = sVar;
}

void print_debug_string(Cstring& s)
{
    double      dTemp=-1.0;
    if( s_sDebugEnvVarName != "" && !bGetEnv(s_sDebugEnvVarName, dTemp) )
        return;

    if(!s_fpDebugFile)
        return;
    
    fprintf(s_fpDebugFile, "%s\n", s.string());
    fflush(s_fpDebugFile);
}

void print_debug_string(char* pc)
{
    double      dTemp=-1.0;
    if( s_sDebugEnvVarName != "" && !bGetEnv(s_sDebugEnvVarName, dTemp) )
        return;

    if(!s_fpDebugFile)
        return;
    
    fprintf(s_fpDebugFile, "%s\n", pc);
    fflush(s_fpDebugFile);
}

void create_debug_file()
{
    double      dTemp=-1.0;
    if( s_sDebugEnvVarName != "" && !bGetEnv(s_sDebugEnvVarName, dTemp) )
        return;

    if( s_fpDebugFile )
        return;
    
    if( s_sDebugFilePath == "" )
        return;

    Cstring     strFileName(s_sDebugFilePath);

    if( bFileExists(strFileName) ) // probably from a previous crash
    {
        nFileAppendVersion(strFileName, false);
    }
    
    s_fpDebugFile = fopen(strFileName.string(), "wt");
}

void close_debug_file()
{
    double      dTemp=-1.0;
    if( s_sDebugEnvVarName != "" && !bGetEnv(s_sDebugEnvVarName, dTemp) )
        return;

    if( !s_fpDebugFile )
        return;
    
    fclose(s_fpDebugFile);
    
    s_fpDebugFile = NULL;
}

void delete_debug_file()
{
    double      dTemp=-1.0;
    if( s_sDebugEnvVarName != "" && !bGetEnv(s_sDebugEnvVarName, dTemp) )
        return;
    
    close_debug_file();

    if( s_fpDebugFile )
        return;  // file is not closed for some reason

    nFileDelete(s_sDebugFilePath);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

// The following bit of code from Bradley Smith is designed to supply the
// IP address of the current machine and work on many different versions of
// Linux

#ifdef __FreeBSD__
#define HAVE_STRUCT_SOCKADDR_SA_LEN 1
#endif

#define MAX_IFS 20

#ifdef LINUX
struct ifreq* 
DTREK_next_ifreq(struct ifreq* req)
{
#ifdef HAVE_STRUCT_SOCKADDR_SA_LEN
  return (struct ifreq*) ((char *)&req->ifr_addr + req->ifr_addr.sa_len);
#else
  return req + 1;
#endif
}
#endif

int
nDTREKGetIPAddress(Cstring *psIPAddress, int a4nIPField[4])
{
  // Return the IP address of this computer as a string "128.1.2.3" and
  // already parsed into an array of integers {128, 1, 2, 3}
  // Return status is 0 on success, and non-zero on failure.
  // On failure, return "128.0.0.1" and {128, 0, 0, 1}

  // From modified algorithm of Bradley Smith

  // TODO: What if there are multiple IP addresses for this computer?
  // TODO: Must be tested on Mac/OSX, Solaris and Windows

  int nStat = -1;
  int i;

#ifdef LINUX

  *psIPAddress = "128.0.0.1"; // Start with a dummy address

  Cstring sTempFile;

  // Try to get IP address of the running computer

  //int get_first_ethernet(int socket_fd, char** name, char** addr) {

  int nSocketFd = -1;
  char* pcName = NULL;
  char* pcAddr = NULL;

  nSocketFd = socket(AF_INET, SOCK_STREAM, 0);

  struct ifreq* if_iter;
  struct ifreq* if_end;
  struct ifconf ifc;
  struct ifreq ifs[MAX_IFS];

  ifc.ifc_len = sizeof(ifs);
  ifc.ifc_req = ifs;

  if (ioctl(nSocketFd, SIOCGIFCONF, &ifc)) {
    nStat = -2;
  } else {
    if_end = (struct ifreq*) (ifc.ifc_buf + ifc.ifc_len);

    for (if_iter = ifc.ifc_req; if_iter < if_end; if_iter = DTREK_next_ifreq(if_iter)) {
      struct ifreq ifreq;

      if (if_iter->ifr_addr.sa_family != AF_INET) {
      continue;
    }

      strncpy(ifreq.ifr_name, if_iter->ifr_name, sizeof(ifreq.ifr_name));
      if (ioctl(nSocketFd, SIOCGIFFLAGS, &ifreq) < 0) {
        continue;
      }
      if (   (ifreq.ifr_flags & (IFF_POINTOPOINT|IFF_LOOPBACK)) == 0
          && (ifreq.ifr_flags & IFF_UP) != 0) {

        if (ioctl(nSocketFd, SIOCGIFADDR, &ifreq)) {
          nStat = -1;
        } else {
          pcName = strdup(ifreq.ifr_name);
          pcAddr = strdup(inet_ntoa(((struct sockaddr_in*)&ifreq.ifr_addr)->sin_addr));

          *psIPAddress = Cstring(pcAddr);

          free (pcName);
          free (pcAddr);
          nStat = 0;
          break;  // Leave the loop
        }
      }
    }
  }

  if (0 != nStat)
    printf("ERROR in nDTREKGetIPAddress: %s (%d)\n", strerror(errno), errno);
  if (0 <= nSocketFd)
    {
      nStat = close(nSocketFd);
      nSocketFd = -1;
      if (0 != nStat)
	printf("ERROR in nDTREKGetIPAddress: %s (%d)\n", strerror(errno), errno);
    }
#endif

#ifdef WIN32

  *psIPAddress = "127.0.0.1"; // Initialize as the standard loopback address

  char* pcName = NULL;
  char* pcAddr = NULL;

  WSADATA wsaData;
  int iRet;
  iRet = WSAStartup(MAKEWORD(1, 1), &wsaData);
  if (iRet != 0) {
	 printf("ERROR in nDTREKGetIPAddress: Error code %d encountered when attempting to initialize Windows Sockets.\n", iRet);
     nStat = -2;
  }
  if (-1 == nStat) {

 	 char ac[80];
	 if (gethostname(ac, sizeof(ac)) == SOCKET_ERROR) {
        printf("ERROR in nDTREKGetIPAddress: Error code %d encountered when getting local host name.\n", WSAGetLastError());
		nStat = -2;
	 }

     if (-1 == nStat) {

        struct hostent *phe = gethostbyname(ac);
		if (phe == 0) {
			printf("ERROR in nDTREKGetIPAddress: Bad host lookup.\n");
			nStat = -2;
		}

		if (-1 == nStat && phe->h_addrtype == AF_INET && phe->h_addr_list[0]) {
			pcName = strdup(phe->h_name);

			struct in_addr addr;
			// IP address will be first element in h_addr_list array
			memcpy(&addr, phe->h_addr_list[0], sizeof(struct in_addr));
			pcAddr = strdup(inet_ntoa(addr));

			*psIPAddress = Cstring(pcAddr);

			free (pcName);
			free (pcAddr);
			nStat = 0;
		}
	}
  }

  WSACleanup();
#endif

  Cstring a4sField[4];

  // There are probably better ways to parse the string

  split(*psIPAddress, a4sField, 4, ".");
  a4nIPField[0] = atoi(a4sField[0].string());
  a4nIPField[1] = atoi(a4sField[1].string());
  a4nIPField[2] = atoi(a4sField[2].string());
  a4nIPField[3] = atoi(a4sField[3].string());

  return (nStat);
}
