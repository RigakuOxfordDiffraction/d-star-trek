//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
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
 
//+Include files

#include <stdlib.h>
#include "Cstring.h"
#include "Cimage_header.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

const Cstring Gs_DTREK_OVERWRITE = "DTREK_OVERWRITE";

//+Function prototypes

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

int main(int   argc,
	 char *argv[])
{
  int            nSaveHeaderOut = 1;  // 0 for no, 1 for yes, 2 for always
  int            nStat = 0;
  bool           bHelpRequested = FALSE;

  int            nNumFiles;         // Number of input files
  Cstring       *psInFile;          // Names of input files
  Cimage_header *poHeaderIn;        // Pointer to input header
  Cstring        sOutFile;

  int            nNumOptions;
  Cstring       *psOption;
  Cstring       *psSubOption;
  Cstring        sTemp;
  int            nCount=0, nPos=0;
  int            nVerbose = 0;
  Cimage_header  oHeaderOut;        // Output header


  vDtrekSetModuleName("dtheadermerge");
  vPrintCopyrightInfo();

  //cout << "\ndtheadermerge: Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl;

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;
  nNumFiles = 0;
  nNumOptions = 0;

  // Count the number of input filenames
  int nNumArg;
  for (nNumArg = 0; nNumArg < argc; nNumArg++)
    {
      if ( ('-' != (char) *argv[nNumArg]) && ('+' != (char) *argv[nNumArg]) )
	{
	  nNumFiles++;
	}
    }

  nNumOptions = nNumFiles;
  if (3 <= nVerbose)
    cout << nNumFiles << " input filenames" << endl;

  psInFile         = new Cstring [nNumFiles];
  psOption         = new Cstring [nNumOptions];

  poHeaderIn     = NULL;
  nNumArg = 0;

  // Now parse input filenames and possible field modifiers
  
  int nFile = 0;

  while (nNumArg < argc)
    {
      if ( ('-' != (char) *argv[nNumArg]) && ('+' != (char) *argv[nNumArg]) )
	{
	  psInFile[nFile] = argv[nNumArg];
	  nFile++;
	}
      else
	{
	  if (0 < nFile)
	    {
	      psOption[nFile-1] = psOption[nFile-1] + argv[nNumArg] + ' ';
	    }
	  else
	    {
	      cout << "\nMajor Error!\n" << endl;
	      exit(1);
	    }
	}
      nNumArg++;
    }

  // Processing special options

  nPos = 0;
  nCount = 0;

  while (nPos != -1) {
    nPos = psOption[0].find(' ', nPos+1);
    nCount++;
  }

  psSubOption = new Cstring[nCount];
  nCount = split(psOption[0], psSubOption, nCount, ' ');

  for (nPos = 0; nPos < nCount; nPos++)
    {
      if ((Cstring) "-w" == psSubOption[nPos])
	{
	  nSaveHeaderOut = 2;
	}
      else if ((Cstring) "-h" == psSubOption[nPos])
	{
	  bHelpRequested = TRUE;
	}
      else if ((Cstring) "-v" == psSubOption[nPos])
	{
	  nVerbose++;
	}
    }

  if (TRUE == bHelpRequested)
    {
      vError(1, "\nHelp requested:\n");
    }
  else if (2 >= argc)
    {
      vError(1, "Not enough arguments!\n");
    }

  delete[] psSubOption;

  cout << "Special Options: " << psOption[0] << endl
       << "Output File: " << psInFile[nNumFiles-1] << endl;

  // Search for first input_file option override

  nPos = 0;
  nCount = 0;

  while (nPos != -1) {
    nPos = psOption[1].find(' ', nPos+1);
    nCount++;
  }

  psSubOption = new Cstring[nCount];
  nCount = split(psOption[1], psSubOption, nCount, ' ');

  if (psSubOption[0] == (Cstring) "-*")
    {
      cout << "Implicit \"+*\" overridden." << endl;
      psOption[1] = psOption[1].after("-*");
    }
  else if (psSubOption[0] == (Cstring) "+*")
    {
      cout << "Note: \"+*\" option unnecessary for first header. "
	   << "(Set by default.)" << endl;
    }
  else
    {
      psOption[1] = "+* " + psOption[1];
    }

  delete[] psSubOption;

  nFile = 1;

  while (nFile < nNumFiles-1)
    {
      // This is the name of an input file (or the name of the outputfile)
      cout << "Reading input file " << psInFile[nFile] << "..." << endl;

      // Check the header input file

      if (bFileExists(psInFile[nFile]))
	{
	  // Read in the header input file
	  
	  poHeaderIn = new Cimage_header(psInFile[nFile]);
	  if (!poHeaderIn->bIsAvailable()) // Opened input file successfully
	    {
	      // Some kind of error with this header, so quit 
	      nSaveHeaderOut = 0;
	      cout << "\nError: Cannot read " << psInFile[nFile]
		   << "; skipping this file..." << endl;
	    }
	  else
	    {
	      // OK, have header
	      // Loop through options and set them in the
	      // header (put them there if they do not exist)
	      // Also make sure they exist in output header.
	      
	      Cimage_header oHeaderTemp;
	      
	      if (3 <= nVerbose)
		{
		  cout << "Opened " << psInFile[nFile] 
		       << " successfully" << endl;
		  cout << "Processing Options...";
		}
		
	      int nPos = 0, nCount = 0;
	      
	      while (nPos != -1) {
		nPos = psOption[nFile].find(' ', nPos+1);
		nCount++;
	      }
	      
	      psSubOption = new Cstring[nCount];
	      
	      nCount = split(psOption[nFile], psSubOption, nCount, ' ');
	      
	      if (3 <= nVerbose)
		{
		  cout << " Found " << nCount << " suboptions: "
		       << psOption[nFile] << endl;
		}
	      
	      for (int i = 0; i < nCount; i++)
		{
		  if ('-' == (char) *psSubOption[i].string())
		    {
		      Cstring sSearchMask = psSubOption[i].after('-');
		      if (3 <= nVerbose)
			{
			  cout << "Keyword mask to delete: \"" 
			       << sSearchMask << "\"" << endl;
			}
		      
		      oHeaderTemp.nDeleteMask(sSearchMask);
		    }
		  else if ('+' == (char) *psSubOption[i].string())
		    {
		      Cstring sSearchMask = psSubOption[i].after('+');
		      if (3 <= nVerbose)
			{
			  cout << "Keyword mask to add: \"" 
			       << sSearchMask << "\"" << endl;
			}
			
		      oHeaderTemp.nCopyMask(*poHeaderIn, sSearchMask);
		    } // else if
		} // for
	      delete[] psSubOption;
	      
	      // Temp to merge!
	      
	      oHeaderOut.nCopyMask(oHeaderTemp, "*");	  
	      
	    } // else

	  delete poHeaderIn;
	  poHeaderIn   = NULL;
	}
      else // Input file does not exist
	{
	  cout << "Error: Input file " << psInFile[nFile] 
	       << " doesn't exist; skipping this file." << endl;
	  nSaveHeaderOut = 0;
	}
      nFile++;
    }

  // Write the file!

  if (0 != nSaveHeaderOut)  // Will try to write the output file
    {
      cout << "\nWriting output file: " << psInFile[nNumFiles-1] 
	   << "..." << flush;

      if (  (bFileExists(psInFile[nNumFiles-1])) && (2 != nSaveHeaderOut) 
	    && ("" != sGetEnv(Gs_DTREK_OVERWRITE)) )
	{
	  cout << "\n\nError: Output file " << psInFile[nNumFiles-1]
	       << " already exists. (Will not overwrite it by default.)" 
	       << endl;
	}
      else // File doesn't exist; can be written
	{
	  if (bFileExists(psInFile[nNumFiles-1]))
	    {
	      cout << "\nWarning: Output file " 
		   << psInFile[nNumFiles-1] 
                   << " already exists!" << endl;
	      if ("" != sGetEnv(Gs_DTREK_OVERWRITE))
		{
		  cout << "File will be OVERWRITTEN!" << endl;
		}
	      else
		{
		  cout << "Existing file will get a version number!" << endl;
		}
	    }

	  oHeaderOut.nWrite(psInFile[nNumFiles-1]);
	  cout << " Done!" << endl;
	}
    }
  else  // Will not write the output file.
    {
      cout << "\nOutput file not written because of previous error(s)."
	   << endl;
      nStat = 1;
    }

  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
    }
  cout << "dtheadermerge - Usage:\n"
       << "dtheadermerge [-h] [-w] input_file1 [input_file_options]\\\n"
       << "              [ [input_file2 [input_file_options] [...] ] \\\n"
       << "              output_file [output_file_options]\n\n"
       << "Command line  options: Description:\n\n"
       << "-w option     Specifying the -w option tells the program to\n"
       << "              overwrite the output_file if it already exists.\n"
       << "              (By default, dtheadermerge will not overwrite the\n"
       << "              output_file if it exists or if there is an error.)\n\n"

       << "-h option     Displays this help text.\n\n"

       << "input_file1   The name of an input image header file.\n"
       << "              The name cannot begin with a + or - character.\n\n"
       << "input_file_options\n"
       << "              Options to specify which keywords or group of keywords\n"
       << "              input_file are included or ignored. By default all \n"
       << "              keywords in the first input_file are included, while\n"
       << "              all keywords in the other input_files are excluded.\n"
       << "              This option takes the form:\n"
       << "                -image_header_keyword  (no whitespace!)\n"
       << "                +image_header_keyword \n"
       << "              where\n"
       << "                -      means exclude keywords which fit selection\n"
       << "                +      means include keywords which fit selection\n"
       << "                image_header_keyword is either a single complete keyword\n"
       << "                         (ex. RAXIS_COMPRESSION_RATIO)\n"
       << "                       or keyword mask ending in an asterisk\n"
       << "                          (ex. RAXIS*).\n"
       << "                       or keyword mask starting in an asterisk\n"
       << "                          (ex. *_GONIO_VALUES).\n"
       << "                       A keyword mask tells the program to check for all the\n"
       << "                       keywords that begin with the letters before/after the\n"
       << "                       asterisk but that may end/start with different letters.\n\n"
       << "                Example 1: to include the CRYSTAL_UNIT_CELL keyword use:\n"
       << "                         +CRYSTAL_UNIT_CELL\n"
       << "                Example 2: to ignore all keywords starting with ROTATION use:\n"
       << "                         -ROTATION*\n"
       << "                Example 3: to ignore all keywords starting with SCAN\n"
       << "                         but keep SCAN_ROTATION use:\n"
       << "                         -SCAN* +SCAN_ROTATION\n"
       << "                         (Note: in this case -SCAN* should come before\n"
       << "                                +SCAN_ROTATION)\n\n"
       << "                INCORRECT: +*SCAN*\n"
       << "                     WRONG! a keyword mask cannot both start and end\n"
       << "                            with an asterisk!\n\n"
       << "   NOTE:  with most shells the * character must be escaped.\n"
       << endl;
  cout << "output_file   The name of the output header file.\n"
       << "              The name cannot begin with a + or - character.  As a\n"
       << "              safety precaution this file will not be written if it\n"
       << "              already exists.\n"
       << "              No command-line options should come after the output_file\n"
       << endl;
  exit (nErrorNum);
}
