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
// DTMainHeaderEdit.cpp  Initial author: RB                27-Jul-2005
//
//  This file contains the member functions of class CDTMainHeaderEdit
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
 
#include "DTMainHeaderEdit.h"

#include <stdio.h>        // for sprintf, until we figure out sstream

#include "Cstring.h"
#include "Cimage.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

int CDTMainHeaderEdit::nExecute(unsigned int argc, char* argv[])
{
    Cstring         sKeyword("");
    Cstring         sValue("");
    Cstring         sLine("");
    Cstring         sTemp("");
    Cstring         sOutput("");

    vDtrekSetModuleName("dtheaderedit");
    vPrintCopyrightInfo();
    //cout << "\n\ndtheaderedit: Copyright (c) 2006 Rigaku\n"
    //     << endl;

    cout << "Enter input image name: ";
    
    getline(cin, sTemp);
    if( !cin.good() || sTemp == "" )
    {
        cout << "dtheaderedit: ABORTING: invalid name." << endl;
        return 1;
    }

    //if( cin.peek() == '\n' )
        //cin.ignore();
    
    cout << "Enter output image name: ";
    
    getline(cin, sOutput);
    if( !cin.good() || sOutput == "" )
    {
        cout << "dtheaderedit: ABORTING: invalid name." << endl;
        return 1;
    }

    //if( cin.peek() == '\n' )
        //cin.ignore();

    Cimage*         poImage = new Cimage (sTemp);

    if (!poImage->bIsAvailable())
    {
        cout << "Error reading input image.\n";
        return 1;
    }

    poImage->m_oHeader.nList();

    sTemp = "Header edited by dtheaderedit";

    poImage->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sTemp);

    int         nStat = 0;
    while( 0 == nStat )
    {
#ifdef WIN32
        cout << "\nEnter 'keyword=value;' pair (or to exit enter Ctrl+Z): ";
#else
        cout << "\nEnter 'keyword=value;' pair (or to exit enter Ctrl+D): ";
#endif
        getline(cin, sLine);
        
        if( cin )
        {
            sValue = sLine.after('=');
            if ("" == sValue)
            {
                // No equal sign found, so try a space charactero
                sValue   = sLine.after(' ');
                sKeyword = sLine.before(' ');
            }
            else
            {
                sKeyword = sLine.before('=');
            }
                sValue   = sValue.before(';');

            if ("" == sKeyword)
            {
                cout << "dtheaderedit: Error reading keyword!" << endl;
                return 2;
            }

            // Remove whitespace at beginning of sKeyword
            while( (' '  == sKeyword.GetAt(0))
                || ('\n' == sKeyword.GetAt(0)) 
                || ('\t' == sKeyword.GetAt(0)) )
                sKeyword = sKeyword.after(0);

            if( sKeyword == Cimage_header::ms_sHeaderBytes
                || sKeyword == Cimage_header::ms_sDim
                || sKeyword == Cimage_header::ms_sSize1
                || sKeyword == Cimage_header::ms_sSize2
                || sKeyword == Cimage_header::ms_sDataType
                || sKeyword == Cimage_header::ms_sCompression
                || sKeyword == Cimage_header::ms_sType ) 
            {
                // We should have some override mechanism here!
                cout << sKeyword << " is a reserved keyword. Header not changed!\n";
            }
            else if ("" == sValue)
            {
                poImage->m_oHeader.nDelete(sKeyword);
                poImage->m_oHeader.nList();
            }
            else
            {
		if (sKeyword == Cimage_header::ms_sSaturatedValue)
		  {
		    // This is necessary since the Cimage::nWrite() method
		    // writes the Cimage::m_fSatValue value to the output image
		    // and not the value found in the header.  So we must
		    // change the internal value as well if it is to stick.
		    
		    float fValue;
		    fValue = atof(sValue.string());
		    poImage->vSetSatValue(fValue);
		  }
                poImage->m_oHeader.nReplaceValue(sKeyword, sValue);
                poImage->m_oHeader.nList();

            }
        }
        else if (cin.eof())
        {
            cout << endl;
            return (poImage->nWrite(sOutput));
        }
        else 
        {
            cout << "dtheaderedit: Error reading value for keyword!" << endl;
            nStat = 1;
            // Need to reset cin to accept blank lines.
        }
    }
    
    return nStat;
}

void CDTMainHeaderEdit::vError(const int nErrorNum, const Cstring& sMessage)
{
}
