//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cspacegroup.cc      Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Cspacegroup.
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

#include <string.h>
#include "Dtrek.h"
#include "Crefln.h"
#include "dtreksys.h"            // sTransSymbol and sGetEnv

#define SPACEGROUP_TABLE_DEFINITION
#include "Cspacegroup.h"         // Class definition and prototypes
                                 // Cspacegroup.h includes others

#include "Cspacegroup_table.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


//+Definitions, constants, and initialization of static member variables

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cspacegroup::Cspacegroup()
{
  (void) nInitValues();
}

Cspacegroup::Cspacegroup (const int nNum)
{
  (void) nInitValues();
  if ( (0 < nNum) && (233 > nNum) )
    {
      m_eThe_State = eSpacegroup_known_state;
      m_nNumber    =   nNum;
      m_eThe_Class = eGetClass();
      m_sName      = sNameFromNumber();
      vSetCentFlagFromName();
      nReadEquivPos();
    }
  else
    {
      m_eThe_State =   eSpacegroup_unknown_state;
      m_eThe_Class = eSpacegroup_unknown_class;
      m_nNumber    =   -1;
    }
}

Cspacegroup::~Cspacegroup()
{
  if (NULL != m_pfRotT)
    {
      delete [] m_pfRotT;
      m_pfRotT = NULL;
    }
  if (NULL != m_pfTrans)
    {
      delete [] m_pfTrans;
      m_pfTrans = NULL;
    }
}

int
Cspacegroup::nInitValues()
{
  m_eThe_State             =  eSpacegroup_unknown_state;
  m_nNumber                = -1;
  m_nNumEquiv              = 0;
  m_nNumLaue               = 0;
  m_sName                  = "unknown";
  m_sEquiv                 = "";
  m_sDescription           = "unknown";
  m_nCentFlag              = 0;
  m_pfRotT                 = NULL;
  m_pfTrans                = NULL;
  return (0);
}

int Cspacegroup::nList(const int nFlag)
{
    printf("\nSpacegroup number: %d", m_nNumber);
    printf("\n             name: %s", m_sName.string());

    printf("\nNum. equiv. posns: %d", m_nNumEquiv);
    if( 0 != m_nCentFlag )
    {
        printf("  (without centering)");
        printf("\n        Centering: %c", cGetCentFlag());
    }
    printf("\n");

    if( 1 > nFlag )
    {
        if( NULL != m_pfRotT && NULL != m_pfTrans )
        {
            float *pfTemp1;
            float *pfTemp2;

            pfTemp1 = m_pfRotT;
            pfTemp2 = m_pfTrans;
            
            for(int i = 1; i <= m_nNumEquiv; i++)
            {
                printf("Equival. position %2d:\n", i);
                
                for(int j = 0; j < 3; j++)
                {
                    printf("     ");
                    printf("%6f, %6f, %6f      %6f\n", pfTemp1[0], 
                                                       pfTemp1[1], 
                                                       pfTemp1[2], 
                                                      *pfTemp2++);
                    pfTemp1 += 3;
                }
                printf("\n");
            }
        }
        else
        {
            printf("Equiv. positions: unavailable\n\n");
        }
    }
    fflush(stdout);


    if( eSpacegroup_unknown_state == m_eThe_State )
        return 1;
    else
        return 0;
}

bool Cspacegroup::bIsExtinct(const int nH, const int nK, const int nL)
{
  // From code found in MADNES which originally came from FILME and
  // Bill Bennett.  This uses no trickiness such as looking
  // at the spacegroup equivalences.  Instead, the extinctions are
  // brute-forced encoded for all 230 spacegroups plus some special ones.
  //
  static int s_nNumErrors = 0;

  switch (m_nNumber)
    {

      // Cases with no extinctions:
    case 1:
    case 2:
    case 3:
    case 6:
    case 10:
    case 16:
    case 25:
    case 47:
    case 75:
    case 81:
    case 83:
    case 89:
    case 99:
    case 111:
    case 115:
    case 123:
    case 143:
    case 147:
    case 149:
    case 150:
    case 156:
    case 157:
    case 162:
    case 164:
    case 168:
    case 174:
    case 175:
    case 177:
    case 183:
    case 187:
    case 189:
    case 191:
    case 195:
    case 200:
    case 207:
    case 215:
    case 221:
    case 232:
    case 233:
    case 234:
    case 235:
    case 237:
      return (FALSE);  // No extinctions
//      break;


    case 4:
    case 11:

      //  -  0k0: k=2n  [2(1) parallel to y (nonstandard: 4(2) or 6(3))]

      if (nH == 0 && nL == 0 && (nK % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 5:
    case 8:
    case 12:
    case 21:
    case 35:
    case 65:

      //  -  hkl: h+k=2n  [C centering]

      if ( ((nH+nK) % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 7:
    case 13:
    case 26:

      //  -  h0l: l=2n  [C perpendicular to y]

      if (nK == 0 && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 9:
    case 15:
    case 36:
    case 63:

      //  -  hkl: h+k=2n  [C centering] and h0l: l=2n  [C perpendicular to y]

      if ( ( ((nH+nK) % 2) != 0) ||
           (nK == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 14:

      //  -  h0l: l=2n  [C perpendicular to y] and 0k0: k=2n [2(1) parallel to y]

      if ( (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && nL == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 17:
    case 77:
    case 84:
    case 93:
    case 173:
    case 176:
    case 182:

      //  -  00l: l=2n  [2(1), 4(2) or 6(3) parallel to z]

      if (nH == 0 && nK == 0 && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 18:
    case 90:
    case 113:

      //  -  h00: h=2n; 0k0: k=2n  [2(1) axes parallel to x and y]

      if (nL == 0) {
        if ( (nH == 0 && (nK % 2) != 0) ||
             (nK == 0 && (nH % 2) != 0)     ) return (TRUE);
      }
      return (FALSE);
//      break;


    case 19:
    case 94:
    case 198:
    case 208:

      //  -  h00: h=2n; 0k0: k=2n; 00l: l=2n  [2(1) axes parallel to x, y & z]

      if (nL == 0) {
        if ( (nH == 0 && (nK % 2) != 0) ||
             (nK == 0 && (nH % 2) != 0)     ) return (TRUE);
      }
      else if (nH == 0 && nK == 0 && (nL % 2) != 0)
        return (TRUE);

      return (FALSE);
//      break;


    case 20:

      //  -  hkl: h+k=2n; 00l: l=2n  [C centering plus 2(1) parallel to z]

      if ( ((nH+nK) % 2) != 0 ||
           (nH == 0 && nK == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 22:
    case 42:
    case 69:
    case 196:
    case 202:
    case 209:
    case 216:
    case 225:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice

      if (((nH+nK) % 2) != 0 ||
          ((nK+nL) % 2) != 0     ) return (TRUE);
      return (FALSE);
//      break;


    case 23:
    case 24:
    case 44:
    case 71:
    case 79:
    case 82:
    case 87:
    case 97:
    case 107:
    case 119:
    case 121:
    case 139:
    case 197:
    case 199:
    case 204:
    case 211:
    case 217:
    case 229:

      //  -  hkl: h+k+l=2n  [I centering]

      if (((nH+nK+nL) % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 27:
    case 49:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl l=2n [C perpendicular to x]

      if ( (nK == 0 || nH == 0) && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 28:

      //  -  h0l: h=2n  [a perpendicular to y]

      if (nK == 0 && (nH % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 29:

      //  -  h0l: h=2n  [a perpendicular to y] and 0kl l=2n [C perpendicular to x]

      if ( (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 30:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl k+l=2n [n perpendicular to x]

     if ( (nK == 0 && (nL % 2) != 0) ||
          (nH == 0 && ((nK+nL) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 31:

      //  -  h0l: h+l=2n  [n perpendicular to y]

      if (nK == 0 && ((nH+nL) % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 32:
    case 55:

      //  -  h0l: h=2n  [a perpendicular to y] and 0kl k=2n [b perpendicular to x]

      if ( (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 33:

      //  -  h0l: h=2n  [a perpendicular to y] and 0kl k+l=2n [n perpendicular to x]

      if ( (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 34:
    case 58:

      //  -  h0l: h+l=2n  [n perpendicular to y] and 0kl k+l=2n [n perpendicular to x]

      if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 37:
    case 66:

      //  -  hkl: h+k=2n  [C centering] and h0l: l=2n  [C perpendicular to y]
      //  -  and 0kl l=2n [C perpendicular to x]

      if ( (((nH+nK) % 2) != 0) ||
           ( (nK == 0 || nH == 0) && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 38:

      //  -  hkl: k+l=2n  [A centering]

      if (((nK+nL) % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 39:

      //  -  hkl: k+l=2n  [A centering] and 0kl: k=2n  [b perpendicular to x]

      if ( (((nK+nL) % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 40:

      //  -  hkl: k+l=2n  [A centering] and h0l: h=2n  [a perpendicular to y]

      if ( (((nK+nL) % 2) != 0) ||
           (nK == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 41:

      //  -  hkl: k+l=2n  [A centering] and 0kl: k=2n  [b perpendicular to x]
      //  -  and h0l: h=2n  [a perpendicular to y]

      if ( (((nK+nL) % 2) != 0) || (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 43:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice
      //     and 0kl k+l=4n [d perp. to x] and h0l h+l=4n [d perp. to y]

      if ( (((nH+nK) % 2) != 0) ||
           (((nK+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 4) != 0) ||
           (nK == 0 && ((nH+nL) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 45:
    case 72:

      //  -  hkl: nH+k+l=2n  [I centering] and 0kl: nK=2n  [b perpendicular to x]
      //  -  and h0l: h=2n  [a perpendicular to y]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 46:

      //  -  hkl: nH+k+l=2n  [I centering] and h0l: h=2n  [a perpendicular to y]

      if ( (((nH+nK+nL) %2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 48:
    case 134:
    case 201:
    case 224:

      //  -  h0l: h+l=2n  [n perpendicular to y] and 0kl k+l=2n [n perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z]

      if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 50:
    case 125:

      //  -  h0l: h=2n  [a perpendicular to y] and 0kl nK=2n [b perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z]

      if ( (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 51:

      //  -  hk0: h=2n  [a perpendicular to z]

      if (nL == 0 && (nH % 2) != 0) return (TRUE);
      return (FALSE);
//      break;

    case 52:

      //  -  h0l: h+l=2n  [n perpendicular to y] and 0kl k+l=2n [n perpendicular to x]
      //  -  and hk0: h=2n  [a perpendicular to z]

     if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
          (nH == 0 && ((nK+nL) % 2) != 0) ||
          (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 53:

      //  -  h0l: h+l=2n  [n perpendicular to y] and hk0: h=2n  [a perpendicular to z]

      if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 54:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl l=2n [C perpendicular to x]
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( ((nK == 0 || nH == 0) && (nL % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 56:
    case 138:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl l=2n [C perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z]

      if ( ((nK == 0 || nH == 0) && (nL % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 57:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl nK=2n [b perpendicular to x]

      if ( (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 59:
    case 85:
    case 129:

      //  -  hk0: h+k=2n  [n perpendicular to z]

      if (nL == 0 && ((nH+nK) % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 60:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl nK=2n [b perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z]

      if ( (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 61:
    case 205:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl nK=2n [b perpendicular to x]
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 62:

      //  -  k+l=2n  [n perpendicular to x] and hk0: h=2n  [a perpendicular to z]

      if ( (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 64:

      //  -  hkl: h+k=2n  [C centering] and h0l: l=2n  [C perpendicular to y]
      //  -  and hk0 h=2n [a perpendicular to z]

      if ( (((nH+nK) % 2) != 0)       ||
           (nK == 0 && (nL % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 67:

      //  -  hkl: h+k=2n  [C centering] and hk0: h=2n  [a perpendicular to z]

      if ( (((nH+nK) % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 68:

      //  -  hkl: h+k=2n  [C centering] and h0l: l=2n  [C perpendicular to y]
      //  -  and 0kl: l=2n [C perpendicular to x]
      //  -  and hk0: h=2n [a perpendicular to z]

      if ( (((nH+nK) % 2) != 0) ||
           ( (nK == 0 || nH == 0) && (nL % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 70:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice
      //     and 0kl: k+l=4n [d perp. to x] and h0l h+l=2n [d perp. to y]
      //     and hk0: nH+nK=4n [d perp. to z]

      if ( (((nH+nK) % 2) != 0) ||
           (((nK+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 4) != 0) ||
           (nK == 0 && ((nH+nL) % 4) != 0) ||
           (nL == 0 && ((nH+nK) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 73:
    case 206:

      //  -  hkl: nH+k+l=2n  [I centering]
      //  -  and h0l: l=2n  [C perpendicular to y]
      //  -  and 0kl: nK=2n [b perpendicular to x]
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 74:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( (((nH+nK+nL) %2) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 76:
    case 78:
    case 91:
    case 95:

      //  -  00l: l=4n  [4(1) or 4(3) parallel to z]

      if (nH == 0 && nK == 0 && (nL % 4) != 0) return (TRUE);
      return (FALSE);
//      break;

    case 80:
    case 98:

      //  -  hkl: h+k+l=2n; 00l: l=4n  [I centering + 4(1) anLong z]

      if ( ((nH+nK+nL) % 2) != 0   ||
           (nH == 0 && nK == 0 && (nL % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 86:

      //  -  hk0: h+k=2n  [n perpendicular to z]
      //  -  and 00l: l=2n  [2(1), 4(2) or 6(3) parallel to z]

      if ( (nL == 0 && ((nH+nK) % 2) != 0) ||
           (nH == 0 && nK == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 88:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and hk0: h=2n  [a perpendicular to z]
      //  -  and 00l: l=4n  [4(1) or 4(3) parallel to z]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nL == 0 && (nH % 2) != 0) ||
           (nH == 0 && nK == 0 && (nL % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 92:
    case 96:

      //  -  00l: l=4n; h00: h=2n; 0k0: k=2n [4(1) or 4(3) along z and
      //     2(1) along x, x and y axes interchangable]

      if (nL == 0) {
        if ( (nH == 0 && (nK % 2) != 0) ||
             (nK == 0 && (nH % 2) != 0)   ) return (TRUE);
      }
      else if (nH == 0 && nK == 0 &&
              (nL % 4) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 100:
    case 117:
    case 127:

      //  -  0kl: nK=2n  [b perpendicular to x] and h0l: h=2n [a perpendicular to y]

      if ( (nH == 0 && (nK % 2) != 0) ||
           (nK == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 101:
    case 116:
    case 132:

      //  -  0kl: l=2n  [C perpendicular to x] and h0l: l=2n  [C perpendicular to y]

      if ( (nH == 0 && (nL % 2) != 0)  ||
           (nK == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 102:
    case 118:
    case 136:

      //  -  0kl: k+l=2n  [n perpendicular to x] and h0l: h+l=2n [n perpendicular to y]

      if ( (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nK == 0 && ((nH+nL) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;

    case 103:
    case 124:

      //  -  0kl: l=2n  [C perpendicular to x] and h0l: l=2n  [C perpendicular to y]
      //  -  and hhl: l=2n

      if ( (nH == 0 && (nK % 2) != 0) ||
           (nK == 0 && (nL % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 104:
    case 128:

      //  -  0kl: k+l=2n  [n perpendicular to x] and h0l: h+l=2n [n perpendicular to y]
      //  -  and hhl: l=2n

      if ( (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 105:
    case 112:
    case 131:

      //  -  hhl: l=2n

      if (nK == nH && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 106:
    case 135:

      //  -  0kl: nK=2n  [b perpendicular to x] and h0l: h=2n [a perpendicular to y]
      //  -  and hhl: l=2n

      if ( (nH == 0 && (nK % 2) != 0) ||
           (nK == 0 && (nH % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 108:
    case 120:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and h0l: l=2n  [C perpendicular to y]
      //  -  and 0kl: l=2n [C perpendicular to x]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 109:
    case 122:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and hhl: 2h+l=4n  and nH(-nH)0: h=2n

      if ( (((nH+nK+nL) %2) != 0)  ||
           (nK == nH && ((nH+nH+nL) % 4) != 0) ||
           (nL == 0 && nK == -nH && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 110:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and 0kl: l=2n [C perpendicular to x]
      //  -  and h0l: l=2n [C perpendicular to y]
      //  -  and hhl: 2h+l=4n  and nH(-nH)0: h=2n

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nH == 0 && (nL % 2) != 0) ||
           (nK == 0 && (nL % 2) != 0) ||
           (nK == nH && ((nH+nH+nL) % 4) != 0) ||
           (nL == 0 && nK == -nH && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 114:

      //  -  hhl: l=2n and h00: h=2n; 0k0: k=2n  [2(1) axes parallel to x and y]

      if ( (nK == nH && (nL % 2) != 0)           ||
           (nL == 0 && nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && nK == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 126:

      //  -  h0l: h+l=2n  [n perpendicular to y] and 0kl k+l=2n [n perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z] and hhl: l=2n

      if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 130:

      //  -  h0l: l=2n  [C perpendicular to y] and 0kl l=2n [C perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z] and hhl: l=2n

      if ( (nK == 0 && (nL % 2) != 0) ||
           (nH == 0 && (nL % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 133:

      //  -  h0l: h=2n  [a perpendicular to y] and 0kl nK=2n [b perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z] and hhl: l=2n

      if ( (nK == 0 && (nH % 2) != 0) ||
           (nH == 0 && (nK % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 137:

      //  -  hk0: h+k=2n  [n perpendicular to z] and hhl: l=2n

      if ( (nL == 0 && ((nH+nK) % 2) != 0) ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 140:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and 0kl: l=2n [C perpendicular to x]
      //  -  and h0l: l=2n [C perpendicular to y]
      //  -  and hk0: h+k=2n  [n perpendicular to z]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nH == 0 && (nL % 2) != 0) ||
           (nK == 0 && (nL % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 141:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and hhl: 2h+l=4n
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( (((nH+nK+nL) %2) != 0)              ||
           (nK == nH && ((nH+nH+nL) % 4) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 142:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and 0kl: l=2n [C perpendicular to x]
      //  -  and h0l: l=2n [C perpendicular to y]
      //  -  and hhl: 2h+l=4n
      //  -  and hk0: h=2n  [a perpendicular to z]

      if ( (((nH+nK+nL) %2) != 0)     ||
           (nH == 0 && (nL % 2) != 0) ||
           (nK == 0 && (nL % 2) != 0) ||
           (nK == nH && ((nH+nH+nL) % 4) != 0) ||
           (nL == 0 && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 144:
    case 145:
    case 151:
    case 152:
    case 153:
    case 154:
    case 171:
    case 172:
    case 180:
    case 181:

      //  -  00l: l=3n  [3(1), 3(2), 6(2) or 6(4) parallel to z]

      if (nH == 0 && nK == 0 && (nL % 3) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 146:
    case 148:
    case 155:
    case 160:
    case 166:

      //  -  hkl: -h+k+l=3n  [Hexagonal axes for rhombohedral cell]

      if (((nK+nL-nH) % 3) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 158:
    case 165:
    case 185:
    case 188:
    case 193:

      //  -  h(-h)l: l=2n

      if (nK == -nH && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;

    case 159:
    case 163:
    case 186:
    case 190:
    case 194:

      //  -  hhl: l=2n

      if (nK == nH && (nL % 2) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 161:
    case 167:

      //  -  hkl: -h+k+l=3n  [nHexagonanL axes for rnHombonHedranL //enLnL]
      //  -  and h(-h)l: l=2n

      if ( (((nK+nL-nH) % 3) != 0) ||
           (nK == -nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 169:
    case 170:
    case 178:
    case 179:

      //  -  00l: l=6n  [6(1) or 6(5) parallel to z]

      if (nH == 0 && nK == 0 && (nL % 6) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 184:
    case 192:

      //  -  h(-h)l: l=2n and hhl: l=2n

      if ( (nK == -nH && (nL % 2) != 0)  ||
           (nK == nH && (nL % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 203:
    case 227:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice
      //     and 0kl k+l=4n and hk0 nH+nK=4n and h0l h+l=4n

      if ( (((nH+nK) % 2) != 0) ||
           (((nK+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 4) != 0) ||
           (nL == 0 && ((nH+nK) % 4) != 0) ||
           (nK == 0 && ((nH+nL) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 210:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n; 00l: l=4n; h00: h=4n; 0k0: k=4n
      //     [F centering (two conditions suffice) + 4(1) axis along z,
      //     all axes interchangable]

      if (((nH+nK) % 2) != 0 ||
          ((nK+nL) % 2) != 0) return (TRUE);
      else if (nL == 0) {
	if ( (nH == 0 && (nK % 4) != 0) ||
             (nK == 0 && (nH % 4) != 0) ) return (TRUE);
      }
      else if (nH == 0 && nK == 0 && (nL % 4) != 0) return (TRUE);

      return (FALSE);
//      break;


    case 212:
    case 213:

      //  -  00l: l=4n; h00: h=4n; 0k0: k=4n  [4(1) or 4(3) along z,
      //                                       all axes interchangable}

      if (nL == 0) {
        if ((nH == 0 && (nK % 4) != 0) ||
           ( nK == 0 && (nH % 4) != 0)     ) return (TRUE);
      }
      else if (nH == 0 && nK == 0 &&
               (nL % 4) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 214:

      //  -  hkl: h+k+l=2n; 00l: l=4n; 0k0: k=4n; h00: h=4n
      //     [I centering + 4(1) along z, all axes interchangable]

      if (((nH+nK+nL) %2) != 0)
          return (TRUE);
      else if (nL == 0) {
	if ((nH == 0 && (nK % 4) != 0) ||
	    (nK == 0 && (nH % 4) != 0)     ) return (TRUE);
      }
      else if (nH == 0 && nK == 0 && (nL % 4) != 0) return (TRUE);
      return (FALSE);
//      break;


    case 218:
    case 223:

      //  -  hhl: l=2n and hknH: k=2n and hkk: h=2n

      if ( (nK == nH && (nL % 2) != 0) ||
           (nH == nL && (nK % 2) != 0) ||
           (nK == nL && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 219:
    case 226:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice
      //  -  and hhl: l=2n and hkh: k=2n and hkk: h=2n

      if ( (((nH+nK) % 2) != 0)        ||
           (((nK+nL) % 2) != 0)        ||
           (nK == nH && (nL % 2) != 0) ||
           (nH == nL && (nK % 2) != 0) ||
	   (nK == nL && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 220:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and hhl: 2h+l=4n and hkh: 2h+k=4n and hkk: h+2k=4n

      if ( (((nH+nK+nL) %2) != 0)              ||
           (nK == nH && ((nH+nH+nL) % 4) != 0) ||
           (nH == nL && ((nH+nH+nK) % 4) != 0) ||
           (nK == nL && ((nK+nK+nH) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;

    case 222:

      //  -  h0l: h+l=2n  [n perpendicular to y] and 0kl k+l=2n [n perpendicular to x]
      //  -  and hk0: h+k=2n  [n perpendicular to z] and
      //  -  hhl: l=2n and hkh: k=2n and hkk: h=2n

      if ( (nK == 0 && ((nH+nL) % 2) != 0) ||
           (nH == 0 && ((nK+nL) % 2) != 0) ||
           (nL == 0 && ((nH+nK) % 2) != 0) ||
	   (nK == nH && (nL % 2) != 0)     ||
           (nH == nL && (nK % 2) != 0)     ||
           (nK == nL && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 228:

      //  -  hkl: h+k=2n, k+l=2n, h+l=2n  [F centering]
      //     Note that any two conditions suffice
      //  -  and hhl: l=2n and hkh: k=2n and hkk: h=2n
      //     and 0kl: k+l=4n [d perp. to x] and h0l h+l=2n [d perp. to y]
      //     and hk0: nH+nK=4n [d perp. to z]

      if ( (((nH+nK) % 2) != 0)            ||
	   (((nK+nL) % 2) != 0)            ||
	   (nK == nH && (nL % 2) != 0)     ||
	   (nH == nL && (nK % 2) != 0)     ||
	   (nK == nL && (nH % 2) != 0)     ||
           (nH == 0 && ((nK+nL) % 4) != 0) ||
           (nK == 0 && ((nH+nL) % 4) != 0) ||
           (nL == 0 && ((nH+nK) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 230:

      //  -  hkl: h+k+l=2n  [I centering]
      //  -  and 0kl: l=2n [C perpendicular to x]
      //  -  and h0l: l=2n [C perpendicular to y]
      //  -  and hk0: h=2n  [a perpendicular to z]
      //  -  and hhl: 2h+l=4n and hkh: 2h+k=4n and hkk: h+2k=4n

      if ( (((nH+nK+nL) %2) != 0)
      || (nH == 0 && (nL % 2) != 0)
      || (nK == 0 && (nL % 2) != 0)
      || (nL == 0 && (nH % 2) != 0)
      || (nK == nH && ((nH+nH+nL) % 4) != 0)
      || (nH == nL && ((nH+nH+nK) % 4) != 0)
      || (nK == nL && ((nK+nK+nH) % 4) != 0) ) return (TRUE);
      return (FALSE);
//      break;


    case 231:

      //  -  hkl: k+l=2n; h00: h=2n  [A2(1)22, nonstandard choice of axes
      //                              for C222(1)]

      if ( ((nK+nL) % 2) != 0 ||
           (nK == 0 && nL == 0 &&
           (nH % 2) != 0             ) ) return (TRUE);
      return (FALSE);
//      break;


    case 236:
    case 238:

      //  -  hhl: l=2n and hkh: k=2n and hkk: h=2n [RnHombonHedranL axes for R3(-)//]

      if ( (nK == nH && (nL % 2) != 0) ||
           (nH == nL && (nK % 2) != 0) ||
           (nK == nL && (nH % 2) != 0) ) return (TRUE);
      return (FALSE);
//      break;

    default:
      s_nNumErrors++;
      if (50 > s_nNumErrors)
	cout << "ERROR in Cspacegroup::bExtinct, unsupported spacegroup number: "
	     << m_nNumber << endl;
#ifdef SSI_PC
      cout << flush;
#endif
      return (FALSE);
//      break;
    }

/*
// --- USE STANDARD SPA//E-GROUP NUMBERS (EXTENDED BEYOND 230 FOR NON-
//     STANDARD SETTINGS)
//
2      GOTO (  1,  1,  1,  4,  5,  1,  7,  5,  9,  1,  4,  5,  7, 14,  9,
//  SG NR:     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
//
      1        1, 17, 18, 19, 20,  5, 22, 23, 23,  1,  7, 27, 28, 29, 30,
//  SG NR:    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
//
      2       31, 32, 33, 34,  5,  9, 37, 38, 39, 40, 41, 22, 43, 23, 45,
//  SG NR:    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
//
      3       46,  1, 48, 27, 50, 51, 52, 53, 54, 32, 56, 57, 34, 59, 60,
//  SG NR:    46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
//
      4       61, 62,  9, 64,  5, 37, 67, 68, 22, 70, 23, 45, 73, 74,  1,
//  SG NR:    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
//
      5       76, 17, 76, 23, 80,  1, 23,  1, 17, 59, 86, 23, 88,  1, 18,
//  SG NR:    76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
//
      6       76, 92, 17, 19, 76, 92, 23, 80,  1,100,101,102,103,104,105,
//  SG NR:    91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,
//
      7      106, 23,108,109,110,  1,105, 18,114,  1,101,100,102, 23,108,
//  SG NR:   106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,
//
      8       23,109,  1,103, 50,126,100,104, 59,130,105,101,133, 48,106,
//  SG NR:   121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,
//
      9      102,137, 56, 23,140,141,142,  1,144,144,146,  1,146,  1,  1,
//  SG NR:   136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,
//
      A      144,144,144,144,146,  1,  1,158,159,146,161,  1,159,  1,158,
//  SG NR:   151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,
//
      B      146,161,  1,169,169,144,144, 17,  1,  1, 17,  1,169,169,144,
//  SG NR:   166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,
//
      C      144, 17,  1,184,158,159,  1,158,  1,159,  1,184,158,159,  1,
//  SG NR:   181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,
//
      D       22, 23, 19, 23,  1, 48, 22,203, 23, 61, 73,  1, 19, 22,210,
//  SG NR:   196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,
//
      E       23,212,212,214,  1, 22, 23,218,219,220,  1,222,218, 48, 22,
//  SG NR:   211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,
//
      F      219,203,228, 23,230,231,  1,  1,  1,  1,236,  1,236), IGROUP
//  SG NR:   226,227,228,229,230,231,232,233,234,235,236,237,238
//
*/

}

eSpacegroup_classes
Cspacegroup::eGetClass(void)
{
  // Return the class of the crystal, where

  if (0 == m_nNumber)
    return (eSpacegroup_unknown_class);
  else if ( (1 == m_nNumber) || (2 == m_nNumber) )
    return (eSpacegroup_triclinic_class);
  else if ( (3 <= m_nNumber) && (15 >= m_nNumber) )
    return (eSpacegroup_monoclinic_class);
  else if ( ((16 <= m_nNumber) && (74 >= m_nNumber)) || (231 == m_nNumber) )
    return (eSpacegroup_orthorhombic_class);
  else if ( (75 <= m_nNumber) && (142 >= m_nNumber) )
    return (eSpacegroup_tetragonal_class);
  else if ( (143 <= m_nNumber) && (167 >= m_nNumber) )
    return (eSpacegroup_trigonal_class);
  else if ( (168 <= m_nNumber) && (194 >= m_nNumber) )
    return (eSpacegroup_hexagonal_class);
  else if ( (232 <= m_nNumber) && (238 >= m_nNumber) )
    return (eSpacegroup_rhombohedral_class);
  else if ( (195 <= m_nNumber) && (230 >= m_nNumber) )
    return (eSpacegroup_cubic_class);
  else
    return (eSpacegroup_unknown_class);
}

char
Cspacegroup::cGetClass(void)
{
  // Return the class of the crystal, where

  if (0 == m_nNumber)
    return ('u');
  else if ( (1 == m_nNumber) || (2 == m_nNumber) )
    return ('a');
  else if ( (3 <= m_nNumber) && (15 >= m_nNumber) )
    return ('m');
  else if ( ((16 <= m_nNumber) && (74 >= m_nNumber)) || (231 == m_nNumber) )
    return ('o');
  else if ( (75 <= m_nNumber) && (142 >= m_nNumber) )
    return ('t');
  else if ( (143 <= m_nNumber) && (167 >= m_nNumber) )
    return ('h');
  else if ( (168 <= m_nNumber) && (194 >= m_nNumber) )
    return ('h');
  else if ( (232 <= m_nNumber) && (238 >= m_nNumber) )
    return ('h');
  else if ( (195 <= m_nNumber) && (230 >= m_nNumber) )
    return ('c');
  else
    return ('u');
}

Cstring
Cspacegroup::sGetClass(void)
{
  // Return the class of the crystal, where

  if (0 == m_nNumber)
    return ("unknown");
  else if ( (1 == m_nNumber) || (2 == m_nNumber) )
    return ("triclinic");
  else if ( (3 <= m_nNumber) && (15 >= m_nNumber) )
    return ("monoclinic");
  else if ( ((16 <= m_nNumber) && (74 >= m_nNumber)) || (231 == m_nNumber) )
    return ("orthorhombic");
  else if ( (75 <= m_nNumber) && (142 >= m_nNumber) )
    return ("tetragonal");
  else if (    (146 == m_nNumber) 
	    || (148 == m_nNumber)
	    || (155 == m_nNumber)
	    || (160 == m_nNumber)
	    || (161 == m_nNumber)
	    || (166 == m_nNumber)
	    || (167 == m_nNumber)
	       )
    return ("rhomb/hexagonal");
  else if ( (143 <= m_nNumber) && (194 >= m_nNumber) )
    // TO DO: Divide into trigonal and hexagonal
    return ("trig/hexagonal");
  else if ( (232 <= m_nNumber) && (238 >= m_nNumber) )
    return ("rhomb/hexagonal");
  else if ( (195 <= m_nNumber) && (230 >= m_nNumber) )
    return ("cubic");
  else
    return ("unknown");
}

int Cspacegroup::nReadEquivPos()
{
  // Read and load equivalent positions for the spacegroup
  // NAMEd by member variable m_sName
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // WARNING:
  // The input file does not contain all the symmetry equivalent positions.
  // It has only the (0,0,0) translation for centered cells.
  // It is missing the inversion symmetry of centrosymmetric cells.
  // For example, spacegroup R-3c (# 167) has 36 symmetric equivalent positions
  // but only 6 are listed in the file.  Inversion multiplies this by 2, then
  // the centering translations multiply that by 3 so:  2 * 3 * 6 = 36.
  //
  // Also these are from International Tables Volume 1 (not Volume A), so
  // in some spacegroups the choice of origin is different.
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  int      i, j;
  int      nStat;
  Cstring  sTemp;
  Cstring  sNameRead;
  int      nNumRead;
  int      nNumEquivRead;
  float    fTemp;

  char*    pcSpacegroupTable;
  char*    pcPointer;
  char*    pcTok;
   

  i = strlen(g_pcSpacegroupTable_1_49)+strlen(g_pcSpacegroupTable_50_99)+strlen(g_pcSpacegroupTable_100_149)
      +strlen(g_pcSpacegroupTable_150_199)+strlen(g_pcSpacegroupTable_200_230)+1;
  pcSpacegroupTable = new char[i];
  strcpy(pcSpacegroupTable,g_pcSpacegroupTable_1_49);
  strcat(pcSpacegroupTable,g_pcSpacegroupTable_50_99);
  strcat(pcSpacegroupTable,g_pcSpacegroupTable_100_149);
  strcat(pcSpacegroupTable,g_pcSpacegroupTable_150_199);
  strcat(pcSpacegroupTable,g_pcSpacegroupTable_200_230);


  /*
  // Test code used to determine that the hardcoded spacegroup information is the same as in DTREK_SPACEGROUP_FILE

  ifstream oIn( "DTREK_SPACEGROUP_FILE");
  char* pc;
  i=0;
  for (pc=strtok(pcSpacegroupTable," \n\t");pc;pc=strtok(NULL," \n\t")) {
      oIn >> sTemp;
      if (!(sTemp==pc)) {
          if ((sTemp.length()==strspn(sTemp.string(),".0123456789")) && (strlen(pc)==strspn(pc,".0123456789"))) {
              float f0,f1;
              sscanf(sTemp.string(),"%f",&f0);
              sscanf(pc,"%f",&f1);
              if (f0 != f1) {
                  i = 1;
              };
          };
      };
  };
  */
      if (NULL != m_pfRotT)
	{
	  delete [] m_pfRotT;
	  m_pfRotT = NULL;
	}
      if (NULL != m_pfTrans)
	{
	  delete [] m_pfTrans;
	  m_pfTrans = NULL;
	}

      m_nNumEquiv = 0;
      nStat = 0;

      pcPointer = pcSpacegroupTable;

      do
	{                // do until we read required info or error during read

      pcTok = strtok(pcPointer,"\t\n "); pcPointer = NULL;
      sTemp = "";
      if (pcTok)
          sTemp = pcTok;
      if (!sTemp.length())
	    {
	      // Error reading the file
	      nStat = 1;
	    }
	  else if ("NO." == sTemp)
	    {
	      // Found line beginning with "NO.", then read spacegroup number,
	      //  spacegroup name and number of equivalent positions

          pcTok = strtok(NULL,"\t\n ");
          if (1!=sscanf(pcTok,"%d",&nNumRead))
              nStat++;
          pcTok = strtok(NULL,"\t\n ");
          sNameRead = pcTok;
          pcTok = strtok(NULL,"\t\n ");
          if (1!=sscanf(pcTok,"%f",&fTemp))
              nStat++;
          nNumEquivRead = (int) fTemp;  // It has a decmal point in the file


          if (!nStat) {
//	  cout << "Read: " << nNumRead << " name: "
//               << sNameRead << " numeq: " << nNumEquivRead << endl;
		  if (sNameRead != m_sName)
		    {
		      // We need to skip all these

		      nNumRead = nNumEquivRead * 12;
		      for (i = 0; (i < nNumRead); i++) 
                  pcTok = strtok(NULL,"\t\n ");
		    }
		  else
		    {
		      // Found the start of our spacegroup info
		      // Allocate extra space for non-primitive cells
		      // A, B, C, F, I, R centered?

		      m_nNumEquiv = nNumEquivRead;
		      m_nNumLaue  = abs(nNumEquivRead);

		      // Since the input file does not contain the inversion
		      // operators for the centrosymmetric space groups, create
		      //  them below.  The flag for a centrosymmetric space
		      // group is a negative number of equivalent positions.

		      if (m_nNumEquiv < 0) m_nNumEquiv = -2 * m_nNumEquiv;

		      m_pfRotT    = new float [m_nNumEquiv * 9];
		      m_pfTrans   = new float [m_nNumEquiv * 3];

		      float *pfTempT, *pfTempR, *pfTemp3, *pfTemp4;
		      pfTempT = m_pfTrans;
		      pfTempR = m_pfRotT;
		      for (i = 0; i < abs(nNumEquivRead); i++)
			{
			  for (j = 0; j < 3; j++)
			    {
                  sscanf(strtok(NULL,"\t\n "),"%f",pfTempT);
                  pfTempT++;
                  sscanf(strtok(NULL,"\t\n "),"%f",pfTempR);
                  sscanf(strtok(NULL,"\t\n "),"%f",pfTempR+1);
                  sscanf(strtok(NULL,"\t\n "),"%f",pfTempR+2);
			      pfTempR += 3;
			    }
			}

		      // Add inversion operators for centrosymmetric space
		      // groups.  Just take what we have and multiply by -1.0,
		      // adjust translations to always be positive by adding
		      // 1.0 to negative ones.

		      pfTemp3 = m_pfTrans;
		      pfTemp4 = m_pfRotT;
		      for (i = abs(nNumEquivRead); i < m_nNumEquiv; i++)
			{
			  for (j = 0; j < 3; j++)
			    {
			      *pfTempT     = -*pfTemp3++;
			      if (*pfTempT < 0.0) *pfTempT = *pfTempT + 1.0f;
			      pfTempT++;
			      *pfTempR++   = -*pfTemp4++;
			      *pfTempR++   = -*pfTemp4++;
			      *pfTempR++   = -*pfTemp4++;
			    }
			}
		      // We are done! Or are we?
		      // Transpose the symmetry matrix
		      pfTemp4 = m_pfRotT;
		      for (i = 0; i < m_nNumEquiv; i++)
			{
			  vTranMatND(3, pfTemp4);
			  pfTemp4 = pfTemp4 + 9;
			}
		    }
		}
	    }
	  else if (!pcTok)
	    {
	      cout << "ERROR in Cspacegroup::nReadEquivPos, EOF.\n";
	      nStat = -1;
	    }
	} while ( (NULL == m_pfRotT) && (0 == nStat) );

    delete[] pcSpacegroupTable;

    return (nStat);

}


/*
Function reads in a user specified space group.
Syntax is the same as for a single entry in the DTREK_SPACEGROUP_FILE.

*/

int Cspacegroup::nReadEquivPos(Cstring& sDef) {

  int		nx,ny;
  int		nStat;
  float		fTemp;
  char*		cp;
  char*		cp_ws="\t\n\r ";
  Cstring	sTemp;
  bool		bUseFriedel;


  sTemp=sDef;
  cp=strtok(sTemp.string(),cp_ws);
  if (!cp) return 1;

  if (NULL != m_pfRotT) { delete [] m_pfRotT;  m_pfRotT = NULL; }
  if (NULL != m_pfTrans) { delete [] m_pfTrans;  m_pfTrans = NULL;}

  m_nNumEquiv = 0;
  nStat = 0;
  if ((!cp) || (strcmp(cp,"NO."))) return 1;
  cp=strtok(NULL,cp_ws);
  // Must read in the number for the spacegroup.  Since spacegroup is non-standard user spacegroup,
  // this number should be greater than 230.
  sscanf(cp,"%d",&m_nNumber);
  cp = strtok(NULL,cp_ws);
  cp = strtok(NULL,cp_ws);
  if ((!cp) || (1 != sscanf(cp,"%f",&fTemp))) return 1;

  m_nNumEquiv=(int) fTemp;
  m_nNumLaue  = abs(m_nNumEquiv);
  if (m_nNumEquiv < 0) {
	  bUseFriedel=TRUE;
	  m_nNumEquiv = -2 * m_nNumEquiv;
  } else bUseFriedel=FALSE;

  m_pfRotT    = new float [m_nNumEquiv * 9];
  m_pfTrans   = new float [m_nNumEquiv * 3];

  float *pfTempT, *pfTempR, *pfTemp3, *pfTemp4;
  pfTempT = m_pfTrans;
  pfTempR = m_pfRotT;
  for (nx= 0; nx< m_nNumLaue; nx++){
	for (ny = 0; ny < 3; ny++)  {
		cp=strtok(NULL,cp_ws);
		if ((!cp) || (1!=sscanf(cp,"%f",pfTempT))) return 1;
		pfTempT++;
		cp=strtok(NULL,cp_ws);
		if ((!cp) || (1!=sscanf(cp,"%f",&pfTempR[0]))) return 1;
		cp=strtok(NULL,cp_ws);
		if ((!cp) || (1!=sscanf(cp,"%f",&pfTempR[1]))) return 1;
		cp=strtok(NULL,cp_ws);
		if ((!cp) || (1!=sscanf(cp,"%f",&pfTempR[2]))) return 1;
		pfTempR += 3;
		}
	}
  if (bUseFriedel) {
      pfTemp3 = m_pfTrans;
      pfTemp4 = m_pfRotT;
      for (nx = m_nNumLaue; nx < 2*m_nNumLaue; nx++) {
		for (ny = 0; ny < 3; ny++) {
			*pfTempT     = -*pfTemp3++;
			if (*pfTempT < 0.0) *pfTempT = *pfTempT + 1.0f;
			pfTempT++;
			*pfTempR++   = -*pfTemp4++;
			*pfTempR++   = -*pfTemp4++;
			*pfTempR++   = -*pfTemp4++;
			}
		}
      // We are done! Or are we?
      // Transpose the symmetry matrix
      pfTemp4 = m_pfRotT;
      for (nx = 0; nx < m_nNumEquiv; nx++)	{
		vTranMatND(3, pfTemp4);
		pfTemp4 = pfTemp4 + 9;
		}
  };
  return 0;
}




Cstring Cspacegroup::sNameFromNumber(void)
{
    static Cstring sTemp;

    if (0 >= m_nNumber)
    {
      sTemp = "unknown";
      return (sTemp);
    } else
      return (sTemp=g_cpSpaceGroupNames[min(230-1,max(0,m_nNumber-1))]);
};

int Cspacegroup::nNumberFromNameCaseInsensitive(void)
{
    int nx;
    Cstring sLow;
    char a32cLow[32];
    sLow = m_sName;
    sLow.downcase();
    for (nx=0;nx<230;nx++) {
      // make everything lowercase for the comparison
      memset(a32cLow, 0, sizeof(a32cLow));  // clear just in case
      strcpy(a32cLow, g_cpSpaceGroupNames[nx]);
      int i = 0;
      while (0 != a32cLow[i]) { a32cLow[i] = tolower(a32cLow[i]); i++; }
      if (!strcmp(a32cLow, sLow.string()))
            break;
    };
    if (nx<230)
        return nx+1;
    else
        return -1;
};

int Cspacegroup::nNumberFromName(void)
{
    int nx;

    for (nx=0;nx<230;nx++) {
        if (!strcmp(g_cpSpaceGroupNames[nx],m_sName.string()))
            break;
    };
    if (nx<230)
        return nx+1;
    else
        return -1;
};
int Cspacegroup::nReduceHKL(Crefln *poRefln, int *pnHKLOut,
			    int *pnFPlusMinus)
{
  // Reduce the reflection's HKL to asymmetric unit HKL
  // Returns 0:     Reflection reduced successfully, is acentric
  //       n>0:     Reflection reduced successfully, is centric with
  //                phases restricted to n / 12 * pi and 180 degrees away
  //        -1:     Reflection is systematically absent;
  //        -2:     Error occurred reducing the hkl probably during packing
  //  pnHKLOut:     the reduced values of the indices
  //  pnFPlusMinus: If this is non-NULL and the input equivalent positions are
  //                valid, then returns 1 if input reflection was an F+, and
  //                                   -1 if input reflection was an F-
  //                It DOES NOT return 0, if the reflection is centric.
  //
  //        T
  // h'  = R  * h,  where R is a direct space rotation matrix derived from the
  // -     =    -              symmetry equivalent position.  T means transpose
  //
  //                                           ___
  // We always have the Friedel operator hkl = hkl
  //
  // For systematic absent test:
  // Look at the translation component of the transformation:
  //  t = h dot T,   phase(h') = phase(h) - 2*pi*t
  //      -     -
  // If an hkl is equal to its transformed hkl, then t must be a whole number
  // or the phase won't be the same unless the reflection is systematically
  // absent (has magnitude=0).
  //
  // From Bill Furey's notes to the CSHL Macromolecular crystallography course:
  //
  // For the centricity test consider:
  // A reflection is centric if symmetry conditions dictate that its phase
  // MUST be one of only two possible values, which differ by 180 deg.  A
  // reflection is centric if, when expanding to all of its symmetry related
  // counterparts, a reflection is generated with indices precisely the negative
  // of the initial indices.

  int i, j;   // Loop counters
  int a3nHKL[3], a3nHKLMinus[3];
  int nStat;
  int nPackMax;
  int nPack1;
  int nPackOrig, nPackMinus;
  int nSymNum;
  int nFlipped;
  int nFlippedOut;
  float *pfTempR, *pfTempT, fTemp;
  int nCentMin;
  int nTemp;
  const int *pnHKLIn;
  bool bAbsent, bCentric;
  bool bFriedel;
  int  nCentric;

  bAbsent  = FALSE;
  bCentric = FALSE;
  bFriedel = FALSE;
  pnHKLIn  = poRefln->pnGetHKL();
  nStat    = -2;                       // Assume error until proven otherwise
  nPackMax = -1;
  nCentMin = 100;
  pfTempR  = m_pfRotT;                    // Stupid pointer tricks

  a3nHKLMinus[0] = -pnHKLIn[0];
  a3nHKLMinus[1] = -pnHKLIn[1];
  a3nHKLMinus[2] = -pnHKLIn[2];

  nPackOrig  = poRefln->nPackHKL();
  nPackMinus = poRefln->nPackHKL(a3nHKLMinus);

  for (i = 0; i < m_nNumEquiv; i++)
    {
      for (j = 0; j < 3; j++)
	{

	  a3nHKL[j] = (int)(pfTempR[0] * pnHKLIn[0]   //            '    T
			  + pfTempR[1] * pnHKLIn[1]   // Eqn above h = (R * h)
			  + pfTempR[2] * pnHKLIn[2]); //           -    =   -
	  pfTempR += 3;

//	  a3nHKL[j] = (int)(*pfTempR++ * pnHKLIn[0]   //            '    T
//			  + *pfTempR++ * pnHKLIn[1]   // Eqn above h = (R * h)
//			  + *pfTempR++ * pnHKLIn[2]); //           -    =   -
	}

//    cout << "a3nHKL: " << a3nHKL[0] << ", " << a3nHKL[1] << ", " << a3nHKL[2] << endl;

      // Choose asymmetric unit with largest packed indices
      // Because of the packing algorithm, this is picking
      // largest H, largest K, largest L.

      nPack1 = poRefln->nPackHKL(a3nHKL);

      // Check for systematic absence

      if (nPackOrig == nPack1)
	{
//	  if (   (pnHKLIn[0] != a3nHKL[0])
//	      || (pnHKLIn[1] != a3nHKL[1])
//	      || (pnHKLIn[2] != a3nHKL[2]) )
//	    cout << "PROBLEM with nPack == nPack1 SYSTEMABSENCE!\n";

	  // h' = h

	  pfTempT = &m_pfTrans[3 * i];     // Point to the corrspding trans part
//          fTemp = *pfTempT++ * pnHKLIn[0]  // Eqn above (check T part)
//                + *pfTempT++ * pnHKLIn[1]
//                + *pfTempT++ * pnHKLIn[2];
	  fTemp = pfTempT[0] * pnHKLIn[0]  // Eqn above (check T part)
                + pfTempT[1] * pnHKLIn[1]
	        + pfTempT[2] * pnHKLIn[2];
	  pfTempT += 3;

	  if ( 0.0 != (fTemp - (int) fTemp))
	    bAbsent = TRUE;
	  else if (0 != m_nCentFlag)
	    {
	      // If space group is not primitive, check special centering trans

	      if ( (1 == m_nCentFlag) || (4 == m_nCentFlag) )
		{
		  // A or F centering
		  fTemp = 0.5f * pnHKLIn[1]  +  0.5f * pnHKLIn[2];
		  if (0.0 != (fTemp - (int) fTemp))
		    bAbsent = TRUE;
		}
	      if ( (2 == m_nCentFlag) || (4 == m_nCentFlag) )
		{
		  // B or F centering
		  fTemp = 0.5f * pnHKLIn[0]  +  0.5f * pnHKLIn[2];
		  if (0.0 != (fTemp - (int) fTemp))
		    bAbsent = TRUE;
		}
	      if ( (3 == m_nCentFlag) || (4 == m_nCentFlag) )
		{
		  // C or F centering
		  fTemp = 0.5f * pnHKLIn[0]  +  0.5f * pnHKLIn[1];
		  if (0.0 != (fTemp - (int) fTemp))
		    bAbsent = TRUE;
		}
	      else if (5 == m_nCentFlag)
		{
		  // I centering
		  fTemp = 0.5f * pnHKLIn[0]  +  0.5f * pnHKLIn[1]
                                            +  0.5f * pnHKLIn[2];
		  if (0.0 != (fTemp - (int) fTemp))
		    bAbsent = TRUE;
		}
	      else if (6 == m_nCentFlag)
		{
		  //  R centering
		  fTemp = (2.0f * pnHKLIn[0]  +  pnHKLIn[1] +  pnHKLIn[2]) / 3.0f;
		  if (0.0 !=  (fTemp - (int) fTemp))    // 2/3, 1/3, 1/3
		    bAbsent = TRUE;
		  else
		    {
		      fTemp =  (pnHKLIn[0] +  2.0f * pnHKLIn[1]
                                           +  2.0f * pnHKLIn[2]) / 3.0f;
		      if (0.0 != (fTemp - (int) fTemp))  // 1/3, 2/3, 2/3
			bAbsent = TRUE;
		    }
		}
	    }
	}  // end of systematic absence if statement

      // Check for centric / acentric

      if (nPackMinus == nPack1)
	{
//	  if (   (pnHKLIn[0] != -a3nHKL[0])
//	      || (pnHKLIn[1] != -a3nHKL[1])
//	      || (pnHKLIn[2] != -a3nHKL[2]) )
//	    cout << "PROBLEM with nPack CENTRIC!\n";

	  // Here if h' indices are the exact negative of h, so must be centric

	  bCentric = TRUE;
	  nCentric = i;
	}

      // Make sure h > 0

      nFlipped = 1;
      if (0 > a3nHKL[0])
	{
	  // -h k l

	  a3nHKL[0] = -a3nHKL[0];
	  a3nHKL[1] = -a3nHKL[1];
	  a3nHKL[2] = -a3nHKL[2];
	  nFlipped  = -1;
	}
      else if (0 == a3nHKL[0])
	{
	  // 0 k l, test k index

	  if (0 > a3nHKL[1])
	    {
	      // 0 -k l

	      a3nHKL[1] = -a3nHKL[1];
	      a3nHKL[2] = -a3nHKL[2];
	      nFlipped  = -1;
	    }
	  else if (0 == a3nHKL[1])
	    {
	      // 0 0 l, test l index
	      if (0 > a3nHKL[2])
		{
		  // 0 0 -l
		  a3nHKL[2] = -a3nHKL[2];
		  nFlipped = -1;
		}
	    }
	}
      if (1 != nFlipped)
	nPack1 = poRefln->nPackHKL(a3nHKL);
      if (nPack1 > nPackMax)
	{
	  nPackMax    = nPack1;
	  nSymNum     = i;
	  nFlippedOut = nFlipped;
	}
    }  // end i loop

  if (bCentric)
    {
      // Figure out phase restriction
      // Phase is restricted to (nCentMin/12 * pi)
      //                    and (nCentMin/12 * pi + pi);

      pfTempT = &m_pfTrans[3 * nCentric];  // Point to corresponding trans vector
//      fTemp = *pfTempT++ * pnHKLIn[0]      // Equation above h dot T
//	    + *pfTempT++ * pnHKLIn[1]
//	    + *pfTempT++ * pnHKLIn[2];
	  fTemp = pfTempT[0] * pnHKLIn[0]  // Eqn above (check T part)
                + pfTempT[1] * pnHKLIn[1]
	        + pfTempT[2] * pnHKLIn[2];
	  pfTempT += 3;

      nTemp = (int) (fTemp * 12.0) % 12; // Phase is in units of pi/12
      if (0 >= nTemp) nTemp = nTemp + 12;// If negative, then add 180 deg
      nCentMin = min(nCentMin,nTemp);    // Keep the minimum, goes 1,2...,12
    }

  nPack1 = poRefln->nUnPackHKL(nPackMax, pnHKLOut);
  if (0 > nPack1)
    nStat = -2;
  else if (bAbsent)
    nStat = -1;            // Systematically absent
  else if (bCentric)
    nStat = nCentMin;      // Centric
  else
    nStat = 0;             // Acentric

  if (NULL != pnFPlusMinus)
    {
      if (0 < nFlippedOut)
	*pnFPlusMinus = +1;
      else
	*pnFPlusMinus = -1;
    }
  return (nStat);
}

void
Cspacegroup::vSet(const Cstring &sNameIn)
{
  // Change so that sNameIn can have any case (upper or lower) of characters
  // and the name will get set to the name in the list of known spacegroups

  Cstring s = sNameIn;

  s.downcase();

  if ("unknown" == s || "none" == s)
    {
      m_sName     = "unknown";
      m_nNumber   = 0;
      m_nCentFlag = 0;
      if (NULL != m_pfRotT)
	{
	  delete [] m_pfRotT;
	  m_pfRotT = NULL;
	}
      if (NULL != m_pfTrans)
	{
	  delete [] m_pfTrans;
	  m_pfTrans = NULL;
	}
      m_nNumEquiv = 0;
    }
  else
    {
      m_sName     = sNameIn;
      m_nNumber   = nNumberFromNameCaseInsensitive();
      m_sName     = sNameFromNumber();
      vSetCentFlagFromName();
      nReadEquivPos();
    }
}

void
Cspacegroup::vSet(const int nNumberIn)
{
  if ( (0 >= nNumberIn) || (231 <= nNumberIn) )
    {
      m_nNumber   = 0;
      m_sName     = "unknown";
      m_nCentFlag = 0;
      if (NULL != m_pfRotT)
	{
	  delete [] m_pfRotT;
	  m_pfRotT = NULL;
	}
      if (NULL != m_pfTrans)
	{
	  delete [] m_pfTrans;
	  m_pfTrans = NULL;
	}
      m_nNumEquiv = 0;
    }
  else
    {
      m_nNumber = nNumberIn;
      m_sName   = sNameFromNumber();
      vSetCentFlagFromName();
      nReadEquivPos();
    }
}

char
Cspacegroup::cGetCentFlag(void)
{
  if (0 == m_nCentFlag)
    return ('P');
  else if (1 == m_nCentFlag)
    return ('C');
  else if (2 == m_nCentFlag)
    return ('C');
  else if (3 == m_nCentFlag)
    return ('C');
  else if (4 == m_nCentFlag)
    return ('F');
  else if (5 == m_nCentFlag)
    return ('I');
  else if (6 == m_nCentFlag)
    return ('R');
  else
    return ('U'); // Unknown
}

Cstring
Cspacegroup::sGetEquivPos(const int nPosn, const int nLaueFlag)
{
  // Return a symmetry equivalent position as a string of form "X,Y,Z".
  // If nLaueFlag != 0, then do not put any translations in the string
  //                    and do not include any inversion (centrosymmetric)
  //                    symmetry.

  int i, j;
  Cstring sTemp;
  sTemp = "";
  float *pfTempR;
  float *pfTempRT;
  float *pfTempT;
  int nNumOps;

  char a3cXYZ[4] = "XYZ";

  nNumOps = m_nNumEquiv;
  if (0 != nLaueFlag)
    nNumOps = m_nNumLaue;
  if ( (0 > nPosn) || (nNumOps <= nPosn) )
    return (sTemp);

  bool bHasSomething;
  pfTempRT = &m_pfRotT[9 * nPosn];   // Matrix is stored as transpose
  pfTempT = &m_pfTrans[3 * nPosn];   //  so access it this way
  for (j = 0; j < 3; j++)
    {
      bHasSomething = FALSE;
      pfTempR = pfTempRT++;
      for (i = 0; i < 3; i++)
	{
	  if (0.0 != *pfTempR)
	    {
	      if (0.0 > *pfTempR)
		{
		  sTemp += '-';
		}
	      else if (bHasSomething)
		{
		  sTemp += '+';
		}
	      bHasSomething = TRUE;
	      sTemp = sTemp + a3cXYZ[i];
	    }
	  pfTempR = pfTempR + 3;
	}
      if (0 == nLaueFlag)
	{
	  // Add translation component if not laue

	  if (-0.001 > *pfTempT)
	    {
	      sTemp += '-';
	    }
	  else if (0.001 < *pfTempT)
	    {
	      sTemp += '+';
	    }
	  if ( (0.49 < fabs(*pfTempT)) && (0.51 > fabs(*pfTempT)) )
	    {
	      sTemp += "1/2";
	    }
	  else if ( (0.24 < fabs(*pfTempT)) && (0.26 > fabs(*pfTempT)) )
	    {
	      sTemp += "1/4";
	    }
	  else if ( (0.74 < fabs(*pfTempT)) && (0.76 > fabs(*pfTempT)) )
	    {
	      sTemp += "3/4";
	    }
	  else if ( (0.332 < fabs(*pfTempT)) && (0.334 > fabs(*pfTempT)) )
	    {
	      sTemp += "1/3";
	    }
	  else if ( (0.666 < fabs(*pfTempT)) && (0.668 > fabs(*pfTempT)) )
	    {
	      sTemp += "2/3";
	    }
	  else if ( (0.166 < fabs(*pfTempT)) && (0.168 > fabs(*pfTempT)) )
	    {
	      sTemp += "1/6";
	    }
	  else if ( (0.833 < fabs(*pfTempT)) && (0.834 > fabs(*pfTempT)) )
	    {
	      sTemp += "5/6";
	    }
	}
      pfTempT++;
      if (j < 2)
	sTemp += ", ";
    }
  return (sTemp);
}

void 
Cspacegroup::vSetCentFlagFromName(void)
{
  m_nCentFlag = 0;
  if ('P' == m_sName.GetAt(0))
    m_nCentFlag = 0;
  else if ('A' == m_sName.GetAt(0))
    m_nCentFlag = 1;
  else if ('B' == m_sName.GetAt(0))
    m_nCentFlag = 2;
  else if ('C' == m_sName.GetAt(0))
    m_nCentFlag = 3;
  else if ('F' == m_sName.GetAt(0))
    m_nCentFlag = 4;
  else if ('I' == m_sName.GetAt(0))
    m_nCentFlag = 5;
  else if ('R' == m_sName.GetAt(0))
    m_nCentFlag = 6;
}
