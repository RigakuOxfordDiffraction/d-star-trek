//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMain.cpp     Initial author:  RB 07-Jan-2005

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
#include "DTMain.h"

#include "Cdetector.h"
#include "Csource.h"
#include "CCommandLine.h"

#if !defined(VC6)
using std::cout;
#endif

const char*     c_pcOut                     = "-out";  
const char*     c_pcRef                     = "-ref";  

static void setargv(char* cmdLine, int* argcPtr, char*** argvPtr);
static char** dupargv(char** argv, int maxargc);
static void freeargv(char **vector);

CDTMain::CDTMain()
{
    vInitialize();
}
///////////////////////////////////////////////////////////////////
CDTMain::CDTMain(unsigned int argc, char* argv[])
{
    vInitialize();

    vSetCommandArguments(argc, argv);
}
//////////////////////////////////////////////////////////////////
CDTMain::~CDTMain()
{
    if( !m_bExternalHeader && m_poHeader )
    {
        delete m_poHeader;
        m_poHeader = NULL;
    }

    vClearCommandArguments();
}
///////////////////////////////////////////////////////////////////
void CDTMain::vInitialize()
{
    m_strCommandLineOptions = "";
    
    m_poHeader = NULL;
    m_bExternalHeader = false;
    
    m_ppcDTREKCommandArguments = NULL;
    m_nDTREKCommandArgumentCount = 0;

    m_sHeaderOut = "";
    m_sReflnListOut = "";
}
///////////////////////////////////////////////////////////////////
void CDTMain::vCreateHeader()
{
    if( m_poHeader )
    {
        printf("WARNING:  Removing old header.\n");
        
        delete m_poHeader;
    }

    m_poHeader = new Cimage_header();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bCreateHeader(Cstring& sFileName)
{
    if( "" == sFileName )
    {
        printf("ERROR: header file name is empty.\n");
        return false;
    }

    if( m_poHeader )
    {
        printf("WARNING:  Removing old header.\n");
        
        delete m_poHeader;
    }

    m_poHeader = new Cimage_header(sFileName);
    
    if( !m_poHeader->bIsAvailable() )
    {
        printf("ERROR: problem loading '%s'\n", sFileName.string());
        
        return false;
    } 
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMain::vSetExternalHeader(Cimage_header* pHeader, bool bDeleteCurrentHeader)
{
    if( bDeleteCurrentHeader && m_poHeader )
    {
        printf("WARNING:  Removing old header.\n");
        
        delete m_poHeader;
    }
    
    m_poHeader = pHeader;
    m_bExternalHeader = true;
}

// Parse -reso option and calculate the resolution range. Return false if error.
bool CDTMain::bParseCalculateResolution(Cstring& sResoCommandLineOption, 
                                        float& fResoMin, 
                                        float& fResoMax,
                                        Cstring& strError)
{
    // Initialize
    fResoMin = -1.0f;
    fResoMax = -1.0f;
    
    strError = "";

    // Parse resolution option parameters
    std::vector<MSCString>      saResoParameters;     

    if( 0 != sResoCommandLineOption.nListToVector(saResoParameters, " ") )
        return false;

    int         nNumberOfParameters = saResoParameters.size();
    
    float      fFirstResolution = -1.0f;
    float      fSecondResolution = -1.0f;

    if( 2 == nNumberOfParameters )  // must be two numbers
    {
        if( saResoParameters[0].bIsNumeric() )
        {
            fFirstResolution = atof(saResoParameters[0].string());
        }
        else
        {
            strError = "-reso option has 2 parameters, but the first one is not numeric.";
            return false;
        }

        if( saResoParameters[1].bIsNumeric() )
        {
            fSecondResolution = atof(saResoParameters[1].string());
        }
        else
        {
            strError = "-reso option has 2 parameters, but the second one is not numeric.";
            return false;
        }

        fResoMin = max(fFirstResolution, fSecondResolution);
        fResoMax = min(fFirstResolution, fSecondResolution);
    
        return true;
    }


    float      fResolutionMaxEdge = -1.0f;
    
    Cstring   sInterpretError("-reso option cannot be interpreted. ");     
    
    if( 1 == nNumberOfParameters ) // must be "edge", "edges", "corner", "corners", or "upperEqualsHighest"
    {
        if( !m_poHeader )
        {
            strError = sInterpretError + Cstring("The image header does not exist.");
            return false;
        }
        
        if( !m_poHeader->bIsAvailable() )
        {
            strError = sInterpretError + Cstring("The image header is not available.");
            return false;
        }

        Cstring    sDetectorResoError(""); 

        if( !bCalculateDetectorResolution(fFirstResolution, fSecondResolution, fResolutionMaxEdge, sDetectorResoError) )
        {
            strError = sInterpretError + sDetectorResoError;
            return false;
        }

        if( saResoParameters[0] == "edge" || saResoParameters[0] == "edges" )
        {
            fResoMin = fFirstResolution;
            fResoMax = fResolutionMaxEdge;
        }
        else if( saResoParameters[0] == "corner" || saResoParameters[0] == "corners" || 
            saResoParameters[0] == "upperEqualsHighest")
        {
            fResoMin = fFirstResolution;
            fResoMax = fSecondResolution;
        }
        else
        {
            strError = "-reso option cannot be interpreted. Unknown parameter:";
            strError += saResoParameters[0];

            return false;
        }

        return true;
    }
    
    strError = "-reso option has wrong syntax. Number of parameters should be 1 or 2.";
    
    return false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse -reso option and calculate the resolution range. Return false if error.
bool CDTMain::bParseCalculateResolution(Cstring& sResoCommandLineOption, 
                                        double& dResoMin, 
                                        double& dResoMax,
                                        Cstring& strError)
{
    // Initialize
    dResoMin = -1.0;
    dResoMax = -1.0;
    
    strError = "";

    // Parse resolution option parameters
    std::vector<MSCString>      saResoParameters;     

    if( 0 != sResoCommandLineOption.nListToVector(saResoParameters, " ") )
        return false;

    int         nNumberOfParameters = saResoParameters.size();
    
    double      dFirstResolution = -1.0f;
    double      dSecondResolution = -1.0f;

    if( 2 == nNumberOfParameters )  // must be two numbers
    {
        if( saResoParameters[0].bIsNumeric() )
        {
            dFirstResolution = atof(saResoParameters[0].string());
        }
        else
        {
            strError = "-reso option has 2 parameters, but the first one is not numeric.";
            return false;
        }

        if( saResoParameters[1].bIsNumeric() )
        {
            dSecondResolution = atof(saResoParameters[1].string());
        }
        else
        {
            strError = "-reso option has 2 parameters, but the second one is not numeric.";
            return false;
        }

        dResoMin = max(dFirstResolution, dSecondResolution);
        dResoMax = min(dFirstResolution, dSecondResolution);
    
        return true;
    }


    double      dResolutionMaxEdge = -1.0;
    
    Cstring   sInterpretError("-reso option cannot be interpreted. ");     
    
    if( 1 == nNumberOfParameters ) // must be "edge", "edges", "corner", "corners", or "upperEqualsHighest"
    {
        if( !m_poHeader )
        {
            strError = sInterpretError + Cstring("The image header does not exist.");
            return false;
        }
        
        if( !m_poHeader->bIsAvailable() )
        {
            strError = sInterpretError + Cstring("The image header is not available.");
            return false;
        }

        Cstring    sDetectorResoError(""); 

        if( !bCalculateDetectorResolution(dFirstResolution, dSecondResolution, dResolutionMaxEdge, sDetectorResoError) )
        {
            strError = sInterpretError + sDetectorResoError;
            return false;
        }

        if( saResoParameters[0] == "edge" || saResoParameters[0] == "edges" )
        {
            dResoMin = dFirstResolution;
            dResoMax = dResolutionMaxEdge;
        }
        else if( saResoParameters[0] == "corner" || saResoParameters[0] == "corners" || 
            saResoParameters[0] == "upperEqualsHighest")
        {
            dResoMin = dFirstResolution;
            dResoMax = dSecondResolution;
        }
        else
        {
            strError = "-reso option cannot be interpreted. Unknown parameter:";
            strError += saResoParameters[0];

            return false;
        }

        return true;
    }
    
    strError = "-reso option has wrong syntax. Number of parameters should be 1 or 2.";
    
    return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void CDTMain::vCommandLineToCommandArgumentArray(char* strCommand)
{
    vClearCommandArguments();
	setargv(strCommand, &m_nDTREKCommandArgumentCount, &m_ppcDTREKCommandArguments);
}

void CDTMain::vCommandLineToCommandArgumentArray(unsigned int argc, char* argv[])
{
    vClearCommandArguments();
	m_nDTREKCommandArgumentCount = argc;
	m_ppcDTREKCommandArguments = dupargv(argv, argc);
}

/////////////////////////////////////////////////////////////////////////////////////////////
void CDTMain::vClearCommandArguments()
{
    if( m_ppcDTREKCommandArguments )
    {
        free((void*)m_ppcDTREKCommandArguments);
        m_ppcDTREKCommandArguments = NULL;
        m_nDTREKCommandArgumentCount = 0;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bCalculateDetectorResolution(float& fFirstResolution, 
                                           float& fSecondResolution, 
                                           float& fResolutionMaxEdge,
                                           Cstring& sError)
{
    // Safety check
    if( !m_poHeader )
        return false; 

    // Initialize
    fFirstResolution = 0.0f;
    fSecondResolution = 0.0f;
    fResolutionMaxEdge = 0.0f;

    sError = "";
    /////////////////////////////////////////////

    Cdetector       oDetector(*m_poHeader);
    if( !oDetector.bIsAvailable() )
    {
        sError = "Detector object is not available.";
        return false;
    }

    Csource         oSource(*m_poHeader);
    if( !oSource.bIsAvailable() )
    {
        sError = "Source object is not available.";
        return false;
    }

    float   a3fS0[3] = {0.0f};
    oSource.vCalcGetS0(a3fS0);            
    
    if( 0 != oDetector.nGetResolution(a3fS0, &fFirstResolution, &fSecondResolution, &fResolutionMaxEdge) )
    {
        sError = "Failed to get resolution from detector object.";
        return false;
    }
    
    double      fWavelength = oSource.fGetWavelength();  
    fFirstResolution    *= fWavelength;   
    fSecondResolution   *= fWavelength;
    fResolutionMaxEdge  *= fWavelength;
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bCalculateDetectorResolution(double& dFirstResolution, 
                                           double& dSecondResolution, 
                                           double& dResolutionMaxEdge,
                                           Cstring& sError)
{
    // Safety check
    if( !m_poHeader )
        return false; 

    // Initialize
    dFirstResolution   = -1.0;
    dSecondResolution  = -1.0;
    dResolutionMaxEdge = -1.0;

    sError = "";
    /////////////////////////////////////////////

    Cdetector       oDetector(*m_poHeader);
    if( !oDetector.bIsAvailable() )
    {
        sError = "Detector object is not available.";
        return false;
    }

    Csource         oSource(*m_poHeader);
    if( !oSource.bIsAvailable() )
    {
        sError = "Source object is not available.";
        return false;
    }

    float   a3fS0[3] = {0.0f};
    oSource.vCalcGetS0(a3fS0);            
    
    //////////////////////////////////////////////////////////////
    float       f1 = -1.0f;
    float       f2 = -1.0f;
    float       f3 = -1.0f;

    if( 0 != oDetector.nGetResolution(a3fS0, &f1, &f2, &f3) )
    {
        sError = "Failed to get resolution from detector object.";
        return false;
    }
    else
    {
        dFirstResolution   = f1;
        dSecondResolution  = f2;
        dResolutionMaxEdge = f3;
    }
    ///////////////////////////////////////////////////////////////

    double      dWavelength = oSource.fGetWavelength();  
    dFirstResolution    *= dWavelength;   
    dSecondResolution   *= dWavelength;
    dResolutionMaxEdge  *= dWavelength;
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copy arguments into class members.
void CDTMain::vSetCommandArguments(unsigned int argc, char* argv[])
{
    vCommandLineToCommandArgumentArray(argc, argv);
}
//////////////////////////////////////////////////////////////////////////////
bool CDTMain::bIsHelpRequest()
{
    if( 1 < m_nDTREKCommandArgumentCount && (0==strcmp(m_ppcDTREKCommandArguments[1],"-h") || 
                                             0==strcmp(m_ppcDTREKCommandArguments[1],"-help")) )

        return true;

    return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bGetCommandOption(const char* ccOption, bool bRemove, int* piStart, int* piEnd)
{
    Cstring     strTemp("");

    return bGetCommandOption(ccOption, strTemp, bRemove, piStart, piEnd);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find a command option ccOption, return its start and end indices (optionally), also return the arguments that come with the option
bool CDTMain::bGetCommandOption(const char* ccOption, Cstring& sArgs, bool bRemove, int* piStart, int* piEnd)
{
    // Initialize

    if( piStart )
        *piStart = -1;
    
    if( piEnd )
        *piEnd = -1;
    
    // Internal variables for storing begin and end of option
    int     iStart = -1;
    int     iEnd   = -1;
    
    ///////////////////////////////////////////////////////////////////
    int      ii=0;
    for(ii=0; ii < m_nDTREKCommandArgumentCount; ii++)
    {
        if( 0 == strcmp(m_ppcDTREKCommandArguments[ii], ccOption) )
            break;
    }
    //////////////////////////////////////////////////////////////////

    if( ii == m_nDTREKCommandArgumentCount ) // i.e. the option is not found
        return false;
    else
    {
        if( piStart )
            *piStart = ii;
        
        iStart = ii;
    }


    // Since the option is found, copy all args
    ii++;
    sArgs = "";
    while( ii < m_nDTREKCommandArgumentCount && !bIsDTCommandLineOptionString(m_ppcDTREKCommandArguments[ii]) )
    {
        if( sArgs != "" )
            sArgs += ' ';

        sArgs += m_ppcDTREKCommandArguments[ii];
        
        ii++;
    }
    
    if( piEnd )
        *piEnd = ii - 1;
    
    iEnd = ii - 1;

    if( bRemove )
    {
        if( !bRemoveCommandArgs(iStart, iEnd) )
            return false;
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bGetCommandOption(const char* ccOption, 
                                double& dArg, 
                                bool bRemove, int* piStart, int* piEnd)
{
    Cstring     sArgs("");

    if( !bGetCommandOption(ccOption, sArgs, bRemove, piStart, piEnd) )
        return false;
    
	if( !sArgs.bIsNumeric() )
		return false;
	
    dArg = atof(sArgs.string());
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bGetCommandOption(const char* ccOption, 
                                std::vector<double>& daArgs, 
                                char* pcListSeparators, 
                                bool bRemove, int* piStart, int* piEnd)
{
    daArgs.clear(); // initilalize
    
    Cstring     sArgs("");

    if( !bGetCommandOption(ccOption, sArgs, bRemove, piStart, piEnd) )
        return false;

    sArgs.nListToVector(daArgs, pcListSeparators);

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bGetCommandOption(const char* ccOption, 
                                std::vector<Cstring>& asArgs, 
                                char* pcListSeparators, 
                                bool bRemove, int* piStart, int* piEnd)
{
    asArgs.clear(); // initilalize
    
    Cstring     sArgs("");

    if( !bGetCommandOption(ccOption, sArgs, bRemove, piStart, piEnd) )
        return false;

    sArgs.nListToVector(asArgs, pcListSeparators);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bRemoveCommandOption(const char* ccOption)
{
    // 2DO
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove command arguments specified by the start and end indices. Return false, if the input indices are invalid. 
bool CDTMain::bRemoveCommandArgs(const int iStart, const int iEnd)
{
    Cstring     sCommand("");
    bool        bIndexRangeFound = false;
    for(int ii=0; ii < m_nDTREKCommandArgumentCount; ii++)
    {
        if( ii < iStart || ii > iEnd )
        {
			const bool addQuotes = !strchr(m_ppcDTREKCommandArguments[ii], '"') &&
				strchr(m_ppcDTREKCommandArguments[ii], ' ');
			if(addQuotes) sCommand += "\"";
            sCommand += m_ppcDTREKCommandArguments[ii];
			if(addQuotes) sCommand += "\"";
            sCommand += ' ';
        }
        else
            bIndexRangeFound = true;
    }

    vCommandLineToCommandArgumentArray(sCommand.string());

    return bIndexRangeFound;
}
//////////////////////////////////////////////////////////////////////////////
void CDTMain::vGetCommandArrayString(Cstring& sCommand)
{
    sCommand = "";
    for(int ii=0; ii < m_nDTREKCommandArgumentCount; ii++)
    {
        sCommand += m_ppcDTREKCommandArguments[ii];
        if( ii != m_nDTREKCommandArgumentCount - 1)
            sCommand += ' ';
    }
}

//////////////////////////////////////////////////////////////////////////////
bool CDTMain::bUpdateCommandOption(const char* ccOption, Cstring& sArgs)
{
    int     iStart = -1;
    int     iEnd   = -1;

    Cstring     sTemp("");
    bool    bOptionExists = bGetCommandOption(ccOption, sTemp, false, &iStart, &iEnd);

    if( bOptionExists )
    {
        if( iStart < iEnd )                       // if it exists and does have argument(s),
            bRemoveCommandArgs(iStart + 1, iEnd); // remove the arguments, leaving the option
    }
    else
        iStart = m_nDTREKCommandArgumentCount;  // add input option after the last argument

    Cstring     sCommand("");
    for(int ii=0; ii <= m_nDTREKCommandArgumentCount; ii++)
    {
        if( ii < m_nDTREKCommandArgumentCount )
        {
			const bool addQuotes = !strchr(m_ppcDTREKCommandArguments[ii], '"') &&
				strchr(m_ppcDTREKCommandArguments[ii], ' ');
			if(addQuotes) sCommand += "\"";
            sCommand += m_ppcDTREKCommandArguments[ii];
			if(addQuotes) sCommand += "\"";
            sCommand += ' ';
        }

        if( ii == iStart )
        {
            if( !bOptionExists )
            {
                sCommand += ccOption;
                sCommand += ' ';
            }

            sCommand += sArgs;
            sCommand += ' ';
        }
    }

    vCommandLineToCommandArgumentArray(sCommand.string());

    return true;
}
//////////////////////////////////////////////////////////////////////////////
void CDTMain::vPrependCommandOption(const char* ccOption)
{
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMain::bOpenSocket(int nPort)
{
    cout << "Attempting to open socket on port " << nPort << "\n";
    CDTCoreSocket* socket = new CDTCoreSocket();
    if (socket->ConnectToServer(nPort, "127.0.0.1"))
    {
        cout << "Succeeded in opening socket.\n";
        m_poSocket = socket;
        return true;
    }
    cout << "Failed to open socket.\n";
    return false;
}
//////////////////////////////////////////////////////////////////////////////
void CDTMain::vSetOutputHeaderPath()
{
    Cstring    sArgs("");     

    // See if the user specified the output header path
    if( !bGetCommandOption(c_pcOut, sArgs) )
        return;

    m_sHeaderOut = sDtrekGetPrefix() + sArgs;
}
//////////////////////////////////////////////////////////////////////////////
void CDTMain::vSetOutputReflnListPath()
{
    Cstring    sArgs("");     

    // See if the user specified output refln list path
    if( !bGetCommandOption(c_pcRef, sArgs) )
        return;

    m_sReflnListOut = sDtrekGetPrefix() + sArgs;
}
//////////////////////////////////////////////////////////////////////////////
/*

NAME

	dupargv -- duplicate an argument vector

SYNOPSIS

	char **dupargv (char **vector, maxargc)

DESCRIPTION

	Duplicate an argument vector.  Simply scans through the
	vector, duplicating each argument until the
	terminating NULL is found.

RETURNS

	Returns a pointer to the argument vector if
	successful. Returns NULL if there is insufficient memory to
	complete building the argument vector.

*/
char** dupargv(char **argv, int maxargc)
{
  int argc;
  char **copy;
  
  if (argv == NULL) return NULL;
  
  // the vector
  for (argc = 0; argv[argc] != NULL && argc < maxargc; argc++);
  copy = (char **) malloc ((argc + 1) * sizeof (char *));
  if (copy == NULL) return NULL;
  
  // the strings
  for (argc = 0; argv[argc] != NULL && argc < maxargc; argc++)
  {
      int len = strlen (argv[argc]);
      copy[argc] = (char *) malloc((len + 1) * sizeof (char *));
      if (copy[argc] == NULL)
	  {
		  freeargv (copy);
		  return NULL;
	  }
      strcpy (copy[argc], argv[argc]);
  }
  copy[argc] = NULL;
  return copy;
}
/*
NAME

	freeargv -- free an argument vector

SYNOPSIS

	void freeargv(char **vector)

DESCRIPTION

	Free an argument vector that was built using setargv.  Simply scans
	through the vector, freeing the memory for each argument until the
	terminating NULL is found, and then frees the vector itself.

RETURNS

	No value.

*/
void freeargv(char **vector)
{
  register char **scan;

  if (vector != NULL)
  {
	  for (scan = vector; *scan != NULL; scan++)
	  {
		  free (*scan);
	  }
      free (vector);
  }
}
/*
NAME

	setargv -- build an argument vector from a string

DESCRIPTION

	Given a pointer to a string, parse the string extracting fields
	separated by whitespace and optionally enclosed within either single
	or double quotes (which are stripped off), and build a vector of
	pointers to copies of the string for each field.  The input string
	remains unchanged.

	All of the memory for the pointer array and copies of the string
	is obtained from malloc.  All of the memory can be returned to the
	system with the single function call freeargv, which takes
	argv as it's argument.

	The memory for the argv array is dynamically expanded as necessary.

RETURNS

	Returns a pointer to the argument vector if successful. Returns NULL
	if the input string pointer is NULL or if there is insufficient
	memory to complete building the argument vector.

NOTES

	In order to provide a working buffer for extracting arguments into,
	with appropriate stripping of quotes and translation of backslash
	sequences, we allocate a working buffer at least as long as the input
	string.  This ensures that we always have enough space in which to
	work, since the extracted arg is never larger than the input string.

	If the input is a null string (as opposed to a NULL pointer), then
	setargv returns an argv that has one arg, a null string.

	Argv is always kept terminated with a NULL arg pointer, so it can
	be passed to freeargv at any time, or returned, as appropriate.
*/

void setargv(char *input, int* pargc, char*** pargv)
{
  char *arg;
  char *copybuf;
  int squote = 0;
  int dquote = 0;
  int bsquote = 0;
  int argc = 0;
  int maxargc = 0;
  char** argv = NULL;
  char **nargv;

  if (input != NULL)
  {
	  copybuf = (char *) alloca (strlen (input) + 1);
	  
	  // Use a do while to always execute the loop once.  Always return an
	  // argv, even for null strings.  See NOTES above, test case below.
      do
	  {
		  // Pick off argv[argc]
		  while (isspace(*input))
		  {
			  input++;
		  }
		  if ((maxargc == 0) || (argc >= (maxargc - 1)))
		  {
			  // argv needs initialization, or expansion
			  if (argv == NULL)
			  {
				  maxargc = 8; //adj. param.
				  nargv = (char **) malloc (maxargc * sizeof (char *));
			  }
			  else
			  {
				  maxargc *= 2;
				  nargv = (char **) realloc (argv, maxargc * sizeof (char *));
			  }
			  if (nargv == NULL)
			  {
				  if (argv != NULL)
				  {
					  freeargv (argv);
					  argv = NULL;
				  }
				  break;
			  }
			  argv = nargv;
			  argv[argc] = NULL;
		  }
		  
		  // Begin scanning arg
		  arg = copybuf;
		  while (*input != '\0')
		  {
			  if (isspace(*input) && !squote && !dquote && !bsquote)
			  {
				  break;
			  }
			  else
			  {
				  if (bsquote)
				  {
					  bsquote = 0;
					  *arg++ = *input;
				  }
#ifndef WIN32
				  else if (*input == '\\')
				  {
					  bsquote = 1;
				  }
#endif
				  else if (squote)
				  {
					  if (*input == '\'')
					  {
						  squote = 0;
					  }
					  else
					  {
						  *arg++ = *input;
					  }
				  }
				  else if (dquote)
				  {
					  if (*input == '"')
					  {
						  dquote = 0;
					  }
					  else
					  {
						  *arg++ = *input;
					  }
				  }
				  else
				  {
					  if (*input == '\'')
					  {
						  squote = 1;
					  }
					  else if (*input == '"')
					  {
						  dquote = 1;
					  }
					  else
					  {
						  *arg++ = *input;
					  }
				  }
				  input++;
			  }
		  }
		  *arg = '\0';
		  argv[argc] = strdup (copybuf);
		  if (argv[argc] == NULL)
		  {
			  freeargv (argv);
			  argv = NULL;
			  break;
		  }
		  argc++;
		  argv[argc] = NULL;
		  while (isspace(*input))
		  {
			  input++;
		  }
	  } while (*input != '\0');
  }
  *pargc = argc;
  *pargv = argv;
}
