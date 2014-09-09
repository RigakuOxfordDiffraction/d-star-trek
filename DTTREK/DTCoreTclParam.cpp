//*****************************************************************
//
// Copyright (c) 2004 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTCoreTclParam.cc   Initial author: RB           04-APR-2004
// This file is a DTREK version of SSI core file  CoreTclParam.cpp

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


// Include files

#ifdef WIN32
#include <tchar.h>
#endif

#include "DTCoreTclParam.h"

#define LPCTSTR (const char*)

#ifndef WIN32

#define _tcsicmp strcasecmp
#define _tcscpy strcpy
#define _tcsnccpy strncpy
#define _tcscat strcat

#endif

///////////////////////////////////////////////////////////////////////////////////////
//             Public functions

// Constructors and destructors
CCoreTclParam::CCoreTclParam ()         // default constructor
{
     tclCommand = _T("");

     paramListHead = NULL;
     paramListLast = NULL;
}

CCoreTclParam::CCoreTclParam(const CString command)
{
     tclCommand = command;

     paramListHead = NULL;
     paramListLast = NULL;
}



//This constructor builds the CCoreTclParam object from an argc, argv list
CCoreTclParam::CCoreTclParam(const int argc, char* argv[])
{
	 int i;

	 if (argc < 1)
     {
          TRACE ("Empty command\n");
          return;
     }

     //tclCommand = (char *)malloc (strlen(argv[0])+1);
     //strcpy (tclCommand, argv[0]);
     tclCommand += argv[0];

     paramListHead = NULL;
     paramListLast = NULL;

     CString   name, value;

	 // Loop through the arguments and see if we need to add any
	 // special characters back in
	 for (i=1; i<argc; i++)
	 {
 		  TclShowSpecialChars(argv[i]);
	 }

     //Loop over the argv list.  Note that the number of name value pairs depends
     //on space around the '=', so i is sometimes updated within this for loop.
     for (i=1; i<argc; i++)
     {
          //Look for the =
          char * pEq = strchr (argv[i], '=');
          if (!pEq) //name = value or name =value
          {
               name = argv[i];
               i++; //The value part will be in the next parameter
               if (i >= argc)
               {
                    TRACE ("Format error in command %s = expected after %s\n", 
                         LPCTSTR(tclCommand), argv[i-1]);
                    return;
               }
               if (argv[i][0] != '=')
               {
                    TRACE ("Format error in command %s at %s\n", LPCTSTR(tclCommand), argv[i]);
                    return;
               }
               if (strlen(argv[i]) == 1) //name = value
               {
                    i++; //The value part will be in the next parameter
                    if (i >= argc)
                    {
                         TRACE ("Format error in command: = expected\n");
                         return;
                    }
                    value = argv[i];
               }
               else //name =value
                    value = &(argv[i][1]);
          }
          else //name= value or name=value case
          {
               if (strlen (pEq) == 1) //name= value
               {
                    *pEq = char (0);
                    name = argv[i];
                    i++;
                    if (i >= argc)
                    {
                         TRACE ("Format error in command: value expected after =\n");
                         return;
                    }
                    value = argv[i];
               }
               else //name=value
               {
                    *pEq = char(0);
                    name = argv[i];
                    pEq++;    //Advance to first character in value;
                    value = pEq;
               }

          }
          SetParamString (name, value);
     }    //Next parameters.
}

CCoreTclParam::~CCoreTclParam() 
{
     //Delete the name value pair list
     NameValue *current, *next;
     next = paramListHead;
     while (next)
     {
          current = next;
          next = current->next;
          delete current;
     }

}

// These methods are used to define a parameter list
bool CCoreTclParam::SetParamString (const CString sName, const CString sValue)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     pNV->value = sValue;
     return true;
}

bool CCoreTclParam::SetParamChar (const CString sName, const char cValue)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     CString sValue = cValue;
     pNV->value = sValue;
     return true;
}

bool CCoreTclParam::SetParamDouble (const CString sName, const double dValue, const int precision)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     //CString sValue;
     //sValue.Format(LPCTSTR("%lf"), dValue);
     //pNV->value = sValue;
     char buffer[256];
     char   format[10];
     sprintf (format, "%s.%dlf", "%", precision);
     sprintf (buffer, format, dValue);
     pNV->value = buffer;
     return true;
}

bool CCoreTclParam::SetParamFloat (const CString sName, const float fValue, const int precision)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     //CString sValue;
     //sValue.Format(LPCTSTR("%f"), fValue);
     //pNV->value = sValue;
     char buffer[256];
     char   format[10];
     sprintf (format, "%s.%df", "%", precision);
     sprintf (buffer, format, fValue);
     pNV->value = buffer;
     return true;
}

bool CCoreTclParam::SetParamInt (const CString sName, const int iValue)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     //CString sValue;
     //sValue.Format(LPCTSTR("%i"), iValue);
     //pNV->value = sValue;
     char buffer[256];
     sprintf (buffer, "%d", iValue);
     pNV->value = buffer;
     return true;
}

bool CCoreTclParam::SetParamShort (const CString sName, const short shValue)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     //CString sValue;
     //sValue.Format(LPCTSTR("%i"), shValue);
     //pNV->value = sValue;
     char buffer[256];
     sprintf (buffer, "%d", shValue);
     pNV->value = buffer;
     return true;
}

bool CCoreTclParam::SetParamAddress (const CString sName, const unsigned int pObj)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
          pNV = NewNameValue (sName);
     if (!pNV)
          return false;

     //CString sValue;
     //sValue.Format(LPCTSTR("%i"), shValue);
     //pNV->value = sValue;
     char buffer[256];
     sprintf (buffer, "%x", pObj);
     pNV->value = buffer;
     return true;
}

bool CCoreTclParam::SetParamBool (const CString sName, const CString sValue)
{
    NameValue* pNV = FindNameValue (sName);
    if (!pNV)
        pNV = NewNameValue (sName);
    if (!pNV)
        return false;

	if (_tcsicmp(sValue,"ON") == 0 ||
		_tcsicmp(sValue,"1")  == 0 ||
		_tcsicmp(sValue,"YES")== 0 )
	    pNV->value = "ON";
	else if (_tcsicmp(sValue,"OFF")== 0 ||
		_tcsicmp(sValue,"0") == 0||
		_tcsicmp(sValue,"NO")== 0 )
	    pNV->value = "OFF";
	else
        pNV->value = sValue;
    return true;
}

bool CCoreTclParam::SetParamBool (const CString sName, const bool bValue)
{
    NameValue* pNV = FindNameValue (sName);
    if (!pNV)
        pNV = NewNameValue (sName);
    if (!pNV)
        return false;

	if (bValue)
	    pNV->value = "ON";
	else
        pNV->value = "OFF";
    return true;
}

// This method is used to generate the command

CString CCoreTclParam::ConstructCommand (void)
{
     int        i;
	 char	    temp[10000];
     CString    command;
     NameValue* pNV = paramListHead;
     bool       brace, firstChar;

     command += tclCommand;
     command += ' ';

     while(pNV)
     {
		  _tcscpy(temp,pNV->name);
		  // Put { } braces around names that contain blanks
          // and don't already have them
          for (i=0, brace=false, firstChar=false; i<(int)strlen(temp); i++)
          {
              if (!firstChar && temp[i] == '{') 
              {
                  brace = false;
                  break;
              }
              else if (!firstChar && temp[i] == ' ') 
                  continue;
              else 
                  firstChar = true;
              if (temp[i] == ' ')
              {
                  brace = true;
                  break;
              }
          }
		  pNV->name = temp;
          if (brace) command += "{";
          command += temp;
          if (brace) command += "}";
          command += " = ";

		  _tcscpy(temp,pNV->value);
          // Put { } braces around values that contain blanks 
          // and don't already have them
          for (i=0, brace=false, firstChar=false; i<(int)strlen(temp); i++)
          {
              if (!firstChar && temp[i] == '{') 
              {
                  brace = false;
                  break;
              }
              else if (!firstChar && temp[i] == ' ') 
                  continue;
              else 
                  firstChar = true;
              if (temp[i] == ' ')
              {
                  brace = true;
                  break;
              }
          }
		  pNV->value = temp;
          if (brace) command += "{";
          command += temp;
          if (brace) command += "}";
          command += ' ';
          pNV = pNV->next;
     }

     return command;
}


//These methods are used to access individual parameters
CString CCoreTclParam::GetParamString (const CString sName, bool& found)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = false;
          return _T("");
     }
     found = true;
     return pNV->value;
}

char CCoreTclParam::GetParamChar (const CString sName, bool& found)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = false;
          return char(0);
     }
     found = true;
     return (char)pNV->value.GetAt(0);
}

double CCoreTclParam::GetParamDouble (const CString sName, bool& found, bool& valid)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = valid = false;
          return 0.0;
     }
     found = true;
     if (IsValidDouble(LPCTSTR(pNV->value)))
     {
          valid = true;
          double dValue;
          sscanf (LPCTSTR(pNV->value),"%lf",&dValue);
          return dValue;
     }
     else
     {
          valid = false;
          return 0.0;
     }
}

float CCoreTclParam::GetParamFloat (const CString sName, bool& found, bool& valid)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = valid = false;
          return 0.0;
     }
     found = true;
     if (IsValidDouble(LPCTSTR(pNV->value)))
     {
          valid = true;
          float fValue;
          sscanf (LPCTSTR(pNV->value),"%f",&fValue);
          return fValue;
     }
     else
     {
          valid = false;
          return 0.0;
     }
}

int CCoreTclParam::GetParamInt (const CString sName, bool& found, bool& valid)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = valid = false;
          return 0;
     }
     found = true;
     if (IsValidInt(LPCTSTR(pNV->value)))
     {
          valid = true;
          int iValue;
          sscanf (LPCTSTR(pNV->value),"%d",&iValue);
          return iValue;
     }
     else
     {
          valid = false;
          return 0;
     }
}

short CCoreTclParam::GetParamShort (const CString sName, bool& found, bool& valid)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = valid = false;
          return 0;
     }
     found = true;
     if (IsValidInt(LPCTSTR(pNV->value)))
     {
          short shValue;
          int iValue;
          sscanf (LPCTSTR(pNV->value),"%d",&iValue);
          if ((iValue > -32768) && (iValue < 32768))
          {
               valid = true;
               shValue = iValue;
               return shValue;
          }
          else
          {
               valid = false;
               return 0;
          }
     }
     else
     {
          valid = false;
          return 0;
     }
}

unsigned int CCoreTclParam::GetParamAddress (const CString sName, bool& found, bool& valid)
{
     NameValue* pNV = FindNameValue (sName);
     if (!pNV)
     {
          found = valid = false;
          return 0;
     }
     found = true;
     if (IsValidAddress(LPCTSTR(pNV->value)))
     {
          valid = true;
          unsigned int iValue;
          sscanf (LPCTSTR(pNV->value),"%x",&iValue);
          return iValue;
     }
     else
     {
          valid = false;
          return 0;
     }
}

bool CCoreTclParam::GetParamBool (const CString sName, bool& found, bool& valid)
{
	bool bValue;
    NameValue* pNV = FindNameValue (sName);

    if (!pNV)
    {
         found = valid = false;
         return 0;
    }
    found = true;

	if (_tcsicmp(LPCTSTR(pNV->value),"ON") == 0)
	  bValue = ON;
	else if (_tcsicmp(LPCTSTR(pNV->value),"OFF") == 0)
	  bValue = OFF;
	else
	{
	  valid = false;
	  return 0;
	}

	valid = true;
	return bValue;
}

int CCoreTclParam::CountParams (void)
{
     NameValue* next = paramListHead;
     int count = 0;
     while (next)
     {
          count++;
          next = next->next;
     }
     return count;
}

CString CCoreTclParam::GetParamName(int index)
{
     NameValue* next = paramListHead;

     //Count the NameValue nodes (we already have index 0 so start at 1)
     for (int i=1; i<=index; i++)
     {
          next = next->next;
          if (!next)
               return _T("");
     }

     return next->name;
}

/////////////////////////////////////////////////////////////////////////////
//                       Private Methods
/////////////////////////////////////////////////////////////////////////////
/*
char * CCoreTclParam::copyString (const CString s)
{
     const char *pString = (const char *)(s);
     char *copy = (char *)malloc (strlen(pString)+1);
     strcpy (copy, pString);
     return copy;
}

*/

void CCoreTclParam::PrintParamList(void)
{
     NameValue* next = paramListHead;
     while (next)
     {
          TRACE ("  Name: %s  Value: %s\n", LPCTSTR(next->name), LPCTSTR(next->value));
          next = next->next;
     }
}

NameValue* CCoreTclParam::FindNameValue (const CString sName)
{
     NameValue* next = paramListHead;
     while (next)
     {
          if (next->name == sName)
               return next;
          next = next->next;
     }
     //Nothing found
     return NULL;
}

NameValue* CCoreTclParam::NewNameValue (const CString sName)
{
     NameValue* pNV = new NameValue;
     if (!pNV)
     {
          printf ("Memory allocation failure for NameValue struct");
          return FALSE;
     }

     pNV->next = NULL;

     //Link the new name value pair node.
     if (paramListLast)
     {
          paramListLast->next = pNV;
          paramListLast = pNV;
     }
     else
          paramListLast = pNV;

     if (!paramListHead)
          paramListHead = pNV;

     pNV->name = sName;
     return pNV;
}

bool CCoreTclParam::IsValidDouble (const char * dString)
{
     int nPeriod = 0;
     int n = strlen(dString);
     if (n <= 0)
          return false;
     if (!isdigit(dString[0]) && (dString[0] != '+') && (dString[0] != '-')
          && (dString[0] != '.'))
          return false;
     if (!isdigit(dString[0]) && (n <= 1))
          return false;
     if (dString[0] == '.')
          nPeriod++;
     for (int i=1; i<n; i++)
     {
          if (dString[i] == '.')
          {
               if (nPeriod != 0)
                    return false;
               nPeriod++;
          }
          else if (!isdigit(dString[i]))
               return false;
     }
     return true;
}

bool CCoreTclParam::IsValidInt (const char * dString)
{
     int n = strlen(dString);
     if (n <= 0)
          return false;
     if (!isdigit(dString[0]) && (dString[0] != '+') && (dString[0] != '-'))
          return false;
     if (!isdigit(dString[0]) && (n <= 1))
          return false;
     for (int i=1; i<n; i++)
     {
          if (!isdigit(dString[i]))
               return false;
     }
     return true;
}

bool CCoreTclParam::IsValidAddress (const char * dString)
{
     int n = strlen(dString);
     if (n <= 0)
          return false;
     if (!IsHexDigit(dString[0]))
          return false;
     if (!IsHexDigit(dString[0]) && (n <= 1))
          return false;
     for (int i=1; i<n; i++)
     {
          if (!IsHexDigit(dString[i]))
               return false;
     }
     return true;
}


bool CCoreTclParam::IsHexDigit (const char x)
{
     return isdigit(x) || x == 'a' || x == 'A' ||
                               x == 'b' || x == 'B' ||
                               x == 'c' || x == 'C' ||
                               x == 'd' || x == 'D' ||
                               x == 'e' || x == 'E' ||
                               x == 'f' || x == 'F';
}


//purpose: Set the name of the command, this is use generaly when the default constructor
//   is called to set the name of the command.
bool CCoreTclParam::SetCommand(const CString CommandName)
{
     tclCommand = CommandName;
     return true;
}

CString CCoreTclParam::GetCommand()
{
     return tclCommand;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Functions to parse out the close brace for the tcl arguments
///////////////////////////////////////////////////////////////////////////////////////////

#define HIDDEN_CLOSE_BRACE 22
#define CLOSE_BRACE 125
#define HIDDEN_OPEN_BRACE 24
#define OPEN_BRACE 123
#define HIDDEN_BACK_SLASH 23
#define BACK_SLASH 92

// Parse off each argument and send it to TclHideChars()
int TclHideSpecialChars( char* Strng )
{
    int	    len, start, end, n;
    char    temp[10000], temp2[10000];
    char    NewString[10000];
    bool    BeginBrace;
    bool    mbcs;

    strcpy(temp,Strng);
    start = end = 0;
    len = strlen(Strng);
    _tcscpy(NewString,"");

    while (end < len)
    {
        BeginBrace = false;
        // If the first character of the string is '{' then skip
        // to the inner loop code so we don't accidentally strip
        // out any spaces enclosed in braces e.g. "{   xyz}"
        if (temp[start] == '{')
            BeginBrace = true;

        // Strip out preceding whitespace
        while (!BeginBrace && (temp[start] == ' ' || temp[start] == char(9)))
            start++;

        end = start+1;

        // If after removing whitespace we encounter a '{' then
        // skip to the inner loop to avoid stripping out any
        // spaces enclosed in braces e.g. "   {   xyz}"
        if (temp[start]  == '{')
            BeginBrace = true;

        // Parse out the next argument by either matching surrounding
        // whitespace or matching open and close braces
		while (BeginBrace || (temp[end] != ' ' && temp[end] != char(9) && end < len))
		{
		    // Check for Shift-JIS mbcs characters
		    if ( ((unsigned char)Strng[end-1] > 0x80 && (unsigned char)Strng[end-1] < 0xA0) ||
    		    ((unsigned char)Strng[end-1] > 0xDF && (unsigned char)Strng[end-1] < 0xFD)  ) 
                mbcs = true;
            else
                mbcs = false;

			if (BeginBrace || (temp[end] == '{'  && !mbcs))
			{
                BeginBrace = false;
                // Inner loop ignores whitespace incrementing end until
                // the braces are matched or we hit the end of the string
                int open_cnt = 1;
                end++;
                while (open_cnt != 0 && end < len)
                {
		            // Check for Shift-JIS mbcs characters
		            if ( ((unsigned char)Strng[end-1] > 0x80 && (unsigned char)Strng[end-1] < 0xA0) ||
    		            ((unsigned char)Strng[end-1] > 0xDF && (unsigned char)Strng[end-1] < 0xFD)  ) 
                        mbcs = true;
                    else
                        mbcs = false;
                    if (temp[end] == '{' && !mbcs)
                        open_cnt++;
                    if (temp[end] == '}' && !mbcs)
                        open_cnt--;
                    if (open_cnt != 0)
                        end++;
                }
            }
            if (end < len) end++;
        }

        // Process substring and append to NewString
        n = end - start;
        _tcsnccpy(temp2,&temp[start],n);
        temp2[n] = '\0';
        TclHideChars(temp2);
        start = end + 1;
        _tcscat(NewString,temp2);
        if (end < len) _tcscat(NewString," ");
	}

	_tcscpy(Strng,NewString);
	return 1;
}


// Hide the special characters in a single argument
int TclHideChars( char* Strng )
{
	unsigned int	i, len;


	len = strlen(Strng);
	for (i=1;i < len; i++)
	{
		// Check for Shift-JIS mbcs characters
		if ( ((unsigned char)Strng[i] > 0x80 && (unsigned char)Strng[i] < 0xA0) ||
    		((unsigned char)Strng[i] > 0xDF && (unsigned char)Strng[i] < 0xFD)  ) 
		{
			i++;
			// Bounds checking
			if( i >= len )		// If true, invalid mbcs character
				break;

			if(Strng[i] == (char)OPEN_BRACE)
				Strng[i] = (char)HIDDEN_OPEN_BRACE;

			if(Strng[i] == (char)CLOSE_BRACE)
				Strng[i] = (char)HIDDEN_CLOSE_BRACE;
		}
		// The backslash is the Tcl escape character.  If it precedes a brace,
		// it must be hidden for the Tcl command to succeed
		if(Strng[i] == (char)BACK_SLASH && i < len-1 && 
		  (Strng[i+1] == (char)OPEN_BRACE || Strng[i+1]  == (char)CLOSE_BRACE))
			Strng[i] = (char)HIDDEN_BACK_SLASH;
	}

	return 1;
}

// Put the special characters back in the string.
int TclShowSpecialChars( char *Strng )
{
	unsigned int i;

	for (i=0; i < strlen(Strng); i++)
	{
		if (Strng[i] == (char)HIDDEN_CLOSE_BRACE)
		{
			Strng[i] = (char)CLOSE_BRACE;
		}
		if (Strng[i] == (char)HIDDEN_BACK_SLASH)
		{
			Strng[i] = (char)BACK_SLASH;
		}
        if (Strng[i] == (char)HIDDEN_OPEN_BRACE)
        {
            Strng[i] = (char)OPEN_BRACE;
        }
	}

	return 1;
}

#undef LPCTSTR

