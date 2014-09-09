//
// Copyright (c) 2007 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cspacegroup_check.cpp   Initial author: RB           01-Aug-2007
//
//  This file contains the member functions of class CDTMainScreener
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

#include "Cspacegroup_check.h"

int CSpaceGroupTableInfo::m_nSpaceGroupFreq[] = 
{
    305, 3986,11, 1957,
    273,1,102,22,
    277,5,239,189,
    144,10450,1930,7,
    9,187,3359,86,
    5,1,7,5,
    2,12,1,1,
    242,3,37,9,
    513,14,2,56,
    6,1,5,14,
    47,8,115,3,
    31,5,4,3,
    1,2,9,23,
    15,13,12,101,
    64,30,23,341,
    1261,548,61,96,
    4,12,2,14,
    3,30,4,27,
    8,5,1,47,
    3,7,12,9,
    7,59,6,3,
    37,48,28,98,
    1,4,3,101,
    2,7,1,44,
    2,1,1,1,
    1,4,1,3,
    1,1,2,1,
    6,9,1,1,
    17,68,1,1,
    2,4,1,2,
    12,22,1,8,
    4,1,2,14,
    19,16,3,1,
    2,2,1,17,
    8,3,17,4,
    11,19,10,21,
    10,40,26,122,
    1,5,1,27,
    1,8,23,1,
    4,3,5,21,
    39,1,13,15,
    17,20,36,1,
    14,16,5,1,
    33,1,1,75,
    2,6,1,4,
    1,6,1,1,
    1,15,1,1,
    1,9,1,7,
    1,9,1,1,
    3,15,1,2,
    1,1,1,3,
    36,5,1,1,
    1,3,1,1,
    1,1,7,1,
    18,6,2,4,
    3,1,5,1,
    22,1,1,4,
    8,1
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define DTREK_DSGTI_GET         0x0001
#define DTREK_DSGTI_DELETE      0x0002
static CSpaceGroupTableInfo*   DoSpacegroupTableInfo(DTREK_WORD wCtrl)
{
    static CSpaceGroupTableInfo*    s_pSpaceGroupTableInfo = NULL;

    if( wCtrl == DTREK_DSGTI_GET )
    {
        if( !s_pSpaceGroupTableInfo )
            s_pSpaceGroupTableInfo = new CSpaceGroupTableInfo();
    }
    else if( wCtrl == DTREK_DSGTI_DELETE && s_pSpaceGroupTableInfo )
    {
        delete s_pSpaceGroupTableInfo;
        s_pSpaceGroupTableInfo = NULL;
        
    }
    
    return s_pSpaceGroupTableInfo;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CSpaceGroupTableInfo::CSpaceGroupTableInfo()
{
    m_sTriclinic_table=
    "0  |	-1.AC1		-1.C1 STOP " 
    "1.	    1		        2 "
    ;

    // The unique b axis monoclinic table from the International Tables (Vol. A)
    m_sMonoclinic_table=
    "100     hkl	0kl		h0l		hk0		h00		0k0		00l	|	2/m.AC1				2/m.C1 STOP "
    "1.      -		-		-		-		-		-		-		3,6				        10 "
    "2.      -		-		-		-		-		k		-		4				        11 "
    "3.      -		-		h		-		h		-		-		7.P1a1(C3)			    13.P12/a1(C3) "
    "4.      -		-		h		-		h		k		-		*				        14.P121/a1(C3) "
    "5.      -		-		l		-		-		-		l		7				        13 "
    "6.      -		-		l		-		-		k		l		*				        14 "
    "7.      -		-		hl		-		h		-		l		7.P1n1(C2)			    13.P12/n1(C2) "
    "8.      -		-		hl		-		h		k		l		*				        14.P121/n1(C2) "
    "9.      hk		k		h		hk		h		k		-		5,8				        12 "
    "10.     hk		k		h,l		hk		h		k		l		9				        15 "
    "11.     kl		kl		l		k		-		k		l		5.A121(C2),8.A1m1(C2)	12.A12/m1(C2) "
    "12.     kl		kl		h,l		k		h		k		l		9.A1n1(C2)			    15.A12/n1(C2) "
    "13.     hkl	kl		hl		hk		h		k		l		5.I121(C3),8.I1m1(C3)	12.I12/m1(C3) "
    "14.     hkl 	kl		h,l		hk		h		k		l		9.I1a1(C3)			    15.I12/a1(C3) "
    ;

    // The unique c axis orthorhombic table for the International Tables (Vol. A)
    m_sOrthorhombic_table=
    "200         hkl	    0kl			h0l			hk0			h00			0k0			00l		|	mmm.AC1					mmm.C1 STOP "
    "1.          -			-			-			-			-			-			-			16,25,						47 "
    "2.          -			-			-			-			-			-			l			17  						* "
    "3.          -			-			-			-			-			k			-			17.P2212					* "
    "4.          -			-			-			-			-			k			l			18.P22121					* "
    "5.          -			-			-			-			h			-			-			17.P2122					* "
    "6.          -			-			-			-			h			-			l			18.P21221					* "
    "7.          -			-			-			-			h			k			-			18							* "
    "8.          -			-			-			-			h			k			l			19							* "
    "9.          -			-			-			h			h			-			-			28.Pm2a,26.P21ma			51 "
    "10.         -			-			-			k			-			k			-			28.P2mb,26.Pm21b			51.Pmmb "
    "11.         -			-			-			hk			h			k			-			31.Pm21n|P21mn				59 "
    "12.         -			-			h			-			h			-			-			28,26.P21am					51.Pmam "
    "13.         -			-			h			h			h			-			-			27.P2aa						49.Pmaa "
    "14.         -			-			h			k			h			k			-			29.P21ab					57.Pmab "
    "15.         -			-			h			hk			h			k			-			30.P2an						53.Pman "
    "16.         -			-			l			-			-			-			l			28.P2cm,26					51.Pmcm "
    "17.         -			-			l			h			h			-			l			29.P21ca					57.Pmca "
    "18.         -			-			l			k			-			k			l			32.P2cb						55.Pmcb "
    "19.         -			-			l			hk			h			k			l			33.P21cn					62.Pmcn "
    "20.         -			-			hl			-			h			-			l			31  						59.Pmmm "
    "21.         -			-			hl			h			h			-			l			30.P2na						53 "
    "22.         -			-			hl			k			h			k			l			33.P21nb					62.Pmnb "
    "23.         -			-			hl			hk			h			k			l			34.P2nn						58.Pmnn "
    "24.         -			k			-			-			-			k			-			28.Pbm2,26.Pb21m			51.Pbmm "
    "25.         -			k			-			h			h			k			-			29.Pb21a					57.Pbma "
    "26.         -			k			-			k			-			k			-			27.Pb2b						49.Pbmb "
    "27.         -			k			-			hk			h			k			-			30.Pb2n						53.Pbmn "
    "28.         -			k			h			-			h			k			-			32							55 "
    "29.         -			k			h			h			h			k			-			*							54.Pbaa "
    "30.         -			k			h			k			h			k			-			*							54.Pbab "
    "31.         -			k			h			hk			h			k			-			*							50 "
    "32.         -			k			l			-			-			k			l			29.Pbc21					57 "
    "33.         -			k			l			h			h			k			l			*							61 "
    "34.         -			k			l			k			-			k			l			*							54.Pbcb "
    "35.         -			k			l			hk			h			k			l			*							60 "
    "36.         -			k			hl			-			h			k			l			33.Pbn21					62.Pbnm "
    "37.         -			k			hl			h			h			k			l			*							60.Pbna "
    "38.         -			k			hl			k			h			k			l			*							56.Pbnb "
    "39.         -			k			hl			hk			h			k			l			*							52.Pbnn "
    "40.         -			l			-			-			-			-			l			28.Pc2m,26.Pcm21			51.Pcmm "
    "41.         -			l			-			h			h			-			l			32.Pc2a						55.Pcma "
    "42.         -			l			-			k			-			k			l			29.Pc21b					57.Pcmb "			
    "43.         -			l			-			hk			h			k			l			33.Pc21n					62.Pcmn "
    "44.         -			l			h			-			h			-			l			29							57.Pcam "
    "45.         -			l			h			h			h			-			l			*							54.Pcaa "
    "46.         -			l			h			k			h			k			l			*							61.Pcab "
    "47.         -			l			h			hk			h			k			l			*							60.Pcan "
    "48.         -			l			l			-			-			-			l			27  						49 "
    "49.         -			l			l			h			h			-			l			*							54 "
    "50.         -			l			l			k			-			k			l			*							54.Pccb "
    "51.         -			l			l			hk			h			k			l			*							56 "
    "52.         -			l			hl			-			h			-			l			30.Pcn2						53.Pcnm "
    "53.         -			l			hl			h			h			-			l			*							50.Pcna "
    "54.         -			l			hl			k			h			k			l			*							60.Pcnb "
    "55.         -			l			hl			hk			h			k			l			*							52.Pcnn "
    "56.         -			kl			-			-			-			k			l			31.Pnm21|Pn21m				59.Pnmm "
    "57.         -			kl			-			h			h			k			l			33.Pn21a					62 "
    "58.         -			kl			-			k			-			k			l			30.Pn2b						53.Pnmb "
    "59.         -			kl			-			hk			h			k			l			34.Pn2n						58.Pnmn "
    "60.         -			kl			h			-			h			k			l			33							62.Pnam "
    "61.         -			kl			h			h			h			k			l			*							56.Pnaa "
    "62.         -			kl			h			k			h			k			l			*							60.Pnab "
    "63.         -			kl			h			hk			h			k			l			*							52.Pnan "
    "64.         -			kl			l			-			-			k			l			30							53.Pncm "
    "65.         -			kl			l			h			h			k			l			*							60.Pnca "
    "66.         -			kl			l			k			-			k			l			*							50.Pncb "
    "67.         -			kl			l			hk			h			k			l			*							52.Pncn "
    "68.         -			kl			hl			-			h			k			l			34							58 "
    "69.         -			kl			hl			h			h			k			l			*							52 "
    "70.         -			kl			hl			k			h			k			l			*							52.Pnnb "
    "71.         -			kl			hl			hk			h			k			l			*							48 "
    "72.         hk			k			h			hk			h			k			-			21,35,38.Cm2m|C2mm			65 "
    "73.         hk			k			h			hk			h			k			l			20							* "
    "74.         hk			k			h			h,k			h			k			-			39.Cm2a|C2mb				67 "
    "75.         hk			k			h,l			hk			h			k			l			36,40.C2cm					63 "
    "76.         hk			k			h,l			h,k			h			k			l			41.C2cb						64 "
    "77.         hk			k,l			h			hk			h			k			l			36.Ccm21,40.Cc2m			63.Ccmm "
    "78.         hk			k,l			h			h,k			h			k			l			41.Cc2a						64.Ccmb "
    "79.         hk			k,l			h,l			hk			h			k			l			37							66 "
    "80.         hk			k,l			h,l			h,k			h			k			l			*							68 "
    "81.         hl			l			hl			h			h			-			l			21.B222,38.Bmm2|B2mm,35.Bm2m	65.Bmmm "
    "82.         hl			l			hl			h			h			k			l			20.B2212					* "
    "83.         hl			l			hl			h,k			h			k			l			36.Bm21b,40.B2mb			63.Bmmb "
    "84.         hl			l			h,l			h			h			-			l			39.Bma2|B2cm				67.Bmam|Bmcm "
    "85.         hl			l			h,l			h,k			h			k			l			41.B2cb						64.Bmab "
    "86.         hl			k,l			hl			h			h			k			l			40.Bbm2,36.Bb21m    		63.Bbmm "
    "87.         hl			k,l			hl			h,k			h			k			l			37.Bb2b						66.Bbmb "
    "88.         hl			k,l			h,l			h			h			k			l			41.Bba2						64.Bbcm "
    "89.         hl			k,l			h,l			h,k			h			k			l			*							68.Bbab|Bbcb "
    "90.         kl			kl			l			k			-			k			l			21.A222,38,35.A2mm			65.Ammm "
    "91.         kl			kl			l			k			h			k			l			20.A2122	     			* "
    "92.         kl			kl			l			h,k			h			k			l			40.Am2a,36.A21ma			63.Amma "
    "93.         kl			kl			h,l			k			h			k			l			40,36.A21am					63.Amam "
    "94.         kl			kl			h,l			h,k			h			k			l			37.A2aa						66.Amaa "
    "95.         kl			k,l			l			k			-			k			l			39							67.Abmm|Acmm "
    "96.         kl			k,l			l			h,k			h			k			l			41.Ac2a						64.Abma "
    "97.         kl			k,l			h,l			k			h			k			l			41							64.Acam "
    "98.         kl			k,l			h,l			h,k			h			k			l			*							68.Abaa|Acaa "
    "99.         hkl	    kl			hl			hk			h			k			l			23,24,44					71 "
    "100.        hkl	    kl			hl			h,k			h			k			l			46.Im2a|I2mb				74 "
    "101.        hkl	    kl			h,l			hk			h			k			l			46							74.Imam|Imcm "
    "102.        hkl	    kl			h,l			h,k			h			k			l			45.I2cb						72.Imcb "
    "103.        hkl	    k,l			hl			hk			h			k			l			46.Ibm2|Ic2m				74.Ibmm|Icmm "
    "104.        hkl	    k,l			hl			h,k			h			k			l			45.Ic2a						72.Icma "
    "105.        hkl	    k,l			h,l			hk			h			k			l			45							72 "
    "106.        hkl	    k,l			h,l			h,k			h			k			l			*							73 "
    "107.        hk,hl,kl	kl			hl			hk			h			k			l			22,42						69 "
    "108.        hk,hl,kl	k,l			4hl,h,l		4hk,h,k		4h			4k			4l			43.F2dd						* "
    "109.        hk,hl,kl	4kl,k,l		h,l			4hk,h,k		4h			4k			4l			43.Fd2d						* "
    "110.        hk,hl,kl	4kl,k,l  	4hl,h,l		h,k			4h			4k			4l			43							* "
    "111.        hk,hl,kl	4kl,k,l		4hl,h,l		4hk,h,k		4h			4k			4l			*							70 "
    ;

    m_sTetragonal_table=
    "400         hkl	hk0		0kl		hhl		00l		0k0		hh0	|	4/m.AC1 4/m.C1	4/mmm.AC2		   4/mmm.C2 STOP "
    "1.          -		-		-		-		-		-		-		75,81	83		89,99,111,115   		123 "
    "2.          -		-		-		-		-		k		-		*		*		90,113					* "
    "3.          -		-		-		-		l		-		-		77		84		93						* "
    "4.          -		-		-		-		l		k		-		*		*		94						* "
    "5.          -		-		-		-		4l		-		-		76,78	*		91,95					* "
    "6.          -		-		-		-		4l		k		-		*		*		92,96					* "
    "7.          -		-		-		l		l		-		-		*		*		105,112 				131 "
    "8.          -		-		-		l		l		k		-		*		*		114 					* "
    "9.          -		-		k		-		-		k		-		*		*		100,117					127 "
    "10.         -		-		k		l		l		k		-		*		*		106 					135 "
    "11.         -		-		l		-		l		-		-		*		*		101,116 				132 "
    "12.         -		-		l		l		l		-		-		*		*		103     				124 "
    "13.         -		-		kl		-		l		k		-		*		*		102,118					136 "
    "14.         -		-		kl		l		l		k		-		*		*		104 					128 "
    "15.         -		hk		-		-		-		k		-		*		85		*						129 "
    "16.         -		hk		-		-		l		k		-		*		86		*						* "
    "17.         -		hk		-		l		l		k		-		*		*		*						137 "
    "18.         -		hk		k		-		-		k		-		*		*		*						125 "
    "19.         -		hk		k		l		l		k		-		*		*		*						133 "
    "20.         -		hk		l		-		l		k		-		*		*		*						138 "
    "21.         -		hk		l		l		l		k		-		*		*		*						130 "
    "22.         -		hk		kl		-		l		k		-		*		*		*						134 "
    "23.         -		hk		kl		l		l		k		-		*		*		*						126 "
    "24.         hkl    hk		kl		l		l		k		-		79,82	87		97,107,121,119			139 "
    "25.         hkl    hk		kl		l		4l		k		-		80		*		98						* "
    "26.         hkl    hk		kl		X,l		4l		k		h		*		*		109,122					* "
    "27.         hkl    hk		k,l		l		l		k		-		*		*		108,120					140 "
    "28.         hkl    hk		k,l		X,l		4l		k		h		*		*		110						* "
    "29.         hkl    h,k		kl		l		4l		k		-		*		88		*						* "
    "30.         hkl    h,k		kl		X,l		4l		k		h		*		*		*						141 "
    "31.         hkl    h,k		k,l		X,l		4l		k		h		*		*		*						142 "
    ;

    m_sTrigonal_table=
    "500         hkil	h-h0l	hh-2-hl	000l |  -3.AC1	  -3.C1		   -3m1.AC2	-3m1.C2		   -31m.AC3  -31m.C3 STOP "
    "1.          -		-		-		-		143			147			150,156		164			149,157		162 "
    "2.          -		-		-		3l		144,145		*			152,154		*			151,153		* "
    "3.          -		-		l		l		*			*			*			*			159			163 "
    "4.          -		l		-		l		*			*			158			165			*			* "
    "5.          Y		3hl		3l		3l		146			148			155,160		166			*			* "
    "6.          Y		3hl,l	3l		6l		*			*			161			167			*			* "
    "7.          Z		W		3l		3l		146			148			155,160		166			*			* "
    "8.          Z		W,l		3l		6l		*			*			161			167			*			* "
    ;

    m_sHexagonal_table=
    "600         h-h0l	hh-2-hl	000l |  6/m.AC1	6/m.C1	6/mmm.AC2	   6/mmm.C2 STOP "
    "1.          -		-		-		168,174	175		177,183,189,187		191 "
    "2.          -		-		l		173		176		182					* "
    "3.          -		-		3l		171,172	*		180,181				* "
    "4.          -		-		6l		169,170	*		178,179				* "
    "5.          -		l		l		*		*		186,190				194 "
    "6.          l		-		l		*		*		185,188				193 "
    "7.          l		l		l		*		*		184					192 "
    ;

    m_sCubic_table=
    "700         h0l    hk0		hkl			0kl		hhl		00l	|	m-3.AC1	m-3.C1	m-3m.AC2 m-3m.C2 STOP "
    "1.          -		-		-			-		-		-		195		200		207,215	221 "
    "2.          -		-		-			-		-		l		198		*		208		* "
    "3.          -		-		-			-		-		4l		*		*		213,212	* "
    "4.          -		-		-			-		l		l		*		*		218		223 "
    "5.          l		h		-			k		-		l		*		205		*		* "
    "6.          -		-		-			kl		-		l		*		201		*		224 "
    "7.          -		-		-			kl		l		l		*		*		*		222 "
    "8.          -		-		hkl			kl		l		l		197,199	204		211,217	229 "
    "9.          -		-		hkl			kl		l		4l		*		*		214		* "
    "10.         -		-		hkl			kl		X,l		4l		*		*		220		* "
    "11.         -		-		hkl			k,l		l		l		*		206		*		* "
    "12.         -		-		hkl			k,l		X,l		4l		*		*		*		230 "
    "13.         -		-		hk,hl,kl	k,l		hl		l		196		202		209,216	225 "
    "14.         -		-		hk,hl,kl	k,l		hl		4l		*		*		210		* "
    "15.         -		-		hk,hl,kl	k,l		h,l		l		*		*		219		226 "
    "16.         -		-		hk,hl,kl	4kl,k,l	hl		4l		*		203		*		227 "
    "17.         -		-		hk,hl,kl	4kl,k,l	h,l		4l		*		*		*		228 "
    ;
    
    m_vecLaueClasses.push_back("-1");       // LAUE_1_BAR       = 0
    m_vecLaueClasses.push_back("2/m");      // LAUE_2_OVER_M
    m_vecLaueClasses.push_back("mmm");      // LAUE_MMM
    m_vecLaueClasses.push_back("4/m");      // LAUE_4_OVER_M
    m_vecLaueClasses.push_back("4/mmm");    // LAUE_4_OVER_MMM
    m_vecLaueClasses.push_back("-3");       // LAUE_3_BAR
    m_vecLaueClasses.push_back("-3m1");     // LAUE_3_BAR_M_1
    m_vecLaueClasses.push_back("-31m");     // LAUE_3_BAR_1_M
    m_vecLaueClasses.push_back("6/m");      // LAUE_6_OVER_M
    m_vecLaueClasses.push_back("6/mmm");    // LAUE_6_OVER_MMM
    m_vecLaueClasses.push_back("m-3");      // LAUE_M_3_BAR
    m_vecLaueClasses.push_back("m-3m");     // LAUE_M_3_BAR_M

    m_mapConditions.insert(std::pair<Cstring,enCond>("hkl",enCond_HKL));
    m_mapConditions.insert(std::pair<Cstring,enCond>("0kl",enCond_0KL)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("h0l",enCond_H0L));
    m_mapConditions.insert(std::pair<Cstring,enCond>("hk0",enCond_HK0)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("h00",enCond_H00));
    m_mapConditions.insert(std::pair<Cstring,enCond>("0k0",enCond_0K0)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("00l",enCond_00L)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("hhl",enCond_HHL)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("hh0",enCond_HH0)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("hkil",enCond_HKIL)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("h-h0l",enCond_H_H0L)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("hh-2-hl",enCond_HH_2_HL)); 
    m_mapConditions.insert(std::pair<Cstring,enCond>("000l",enCond_000L));

    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("h",    DTREK_SGC_REFLN_RESTR_H   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("k",    DTREK_SGC_REFLN_RESTR_K   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("l",    DTREK_SGC_REFLN_RESTR_L   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("hk",   DTREK_SGC_REFLN_RESTR_HK  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("kl",   DTREK_SGC_REFLN_RESTR_KL  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("hl",   DTREK_SGC_REFLN_RESTR_HL  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4h",   DTREK_SGC_REFLN_RESTR_4H  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4k",   DTREK_SGC_REFLN_RESTR_4K  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4l",   DTREK_SGC_REFLN_RESTR_4L  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4hk",  DTREK_SGC_REFLN_RESTR_4HK ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4kl",  DTREK_SGC_REFLN_RESTR_4KL ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("4hl",  DTREK_SGC_REFLN_RESTR_4HL ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("hkl",  DTREK_SGC_REFLN_RESTR_HKL ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("W",    DTREK_SGC_REFLN_RESTR_W   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("X",    DTREK_SGC_REFLN_RESTR_X   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("Y",    DTREK_SGC_REFLN_RESTR_Y   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("Z",    DTREK_SGC_REFLN_RESTR_Z   ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("3h",   DTREK_SGC_REFLN_RESTR_3H  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("3l",   DTREK_SGC_REFLN_RESTR_3L  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("6l",   DTREK_SGC_REFLN_RESTR_6L  ));  
    m_mapRestrictions.insert(std::pair<Cstring,DTREK_DWORD>("3hl",  DTREK_SGC_REFLN_RESTR_3HL ));  
                                                                                                         
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("h",  "h=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("k",  "k=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("l",  "l=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("hk", "h+k=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("kl", "k+l=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("hl", "h+l=2n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4h", "h=4n"	  ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4k", "k=4n"	  ));    
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4l", "l=4n"	  ));    
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4hk","h+k=4n"	  ));  
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4kl","k+l=4n"	  ));  
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("4hl","h+l=4n"	  ));  
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("hkl","h+k+l=2n" ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("W",  "l-h=3n"   ));  
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("X",  "2h+l=4n"  )); 
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("Y",  "-h+k+l=3n"));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("Z",  "h-k+l=3n" ));
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("3h", "h=3n"     )); 
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("3l", "l=3n"     )); 
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("6l", "l=6n"     )); 
    m_mapRestrictionLongNames.insert(std::pair<Cstring,Cstring>("3hl","h+l=3n"   ));

    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("h",  "h!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("k",  "k!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("l",  "l!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("hk", "h+k!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("kl", "k+l!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("hl", "h+l!=2n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4h", "h!=4n"	  ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4k", "k!=4n"	  ));    
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4l", "l!=4n"	  ));    
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4hk","h+k!=4n"	  ));  
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4kl","k+l!=4n"	  ));  
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("4hl","h+l!=4n"	  ));  
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("hkl","h+k+l!=2n" ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("W",  "l-h!=3n"   ));  
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("X",  "2h+l!=4n"  )); 
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("Y",  "-h+k+l!=3n"));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("Z",  "h-k+l!=3n" ));
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("3h", "h!=3n"     )); 
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("3l", "l!=3n"     )); 
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("6l", "l!=6n"     )); 
    m_mapRestrictionInvertedLongNames.insert(std::pair<Cstring,Cstring>("3hl","h+l!=3n"   ));

    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_H,   2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_K,   2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_L,   2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_KL,  2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_HL,  2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_HK,  2)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4H,  4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4K,  4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4L,  4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4KL, 4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4HL, 4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_4HK, 4));
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_HKL, 2));
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_W,   3)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_X,   4)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_Y,   3)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_Z,   3)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_3H,  3)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_3L,  3)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_6L,  6)); 
    m_mapRestrictionParity.insert(std::pair<DTREK_DWORD, int>(DTREK_SGC_REFLN_RESTR_3HL, 3));

    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_H,   enRestr_H));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_K,   enRestr_K));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_L,   enRestr_L));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_KL,  enRestr_KL));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_HL,  enRestr_HL));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_HK,  enRestr_HK));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4H,  enRestr_4H));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4K,  enRestr_4K));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4L,  enRestr_4L));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4KL, enRestr_4KL));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4HL, enRestr_4HL));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_4HK, enRestr_4HK));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_HKL, enRestr_HKL));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_W,   enRestr_W));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_X,   enRestr_X));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_Y,   enRestr_Y));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_Z,   enRestr_Z));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_3H,  enRestr_3H));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_3L,  enRestr_3L));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_6L,  enRestr_6L));
    m_mapRestrictionEnum.insert(std::pair<DTREK_DWORD, enRestr>(DTREK_SGC_REFLN_RESTR_3HL, enRestr_3HL));
                                                                                                            
    m_vecRestrictions.push_back(Cstring("h"));
    m_vecRestrictions.push_back(Cstring("k")); 
    m_vecRestrictions.push_back(Cstring("l")); 
    m_vecRestrictions.push_back(Cstring("hk")); 
    m_vecRestrictions.push_back(Cstring("kl"));
    m_vecRestrictions.push_back(Cstring("hl"));
    m_vecRestrictions.push_back(Cstring("4h"));
    m_vecRestrictions.push_back(Cstring("4k"));
    m_vecRestrictions.push_back(Cstring("4l"));
    m_vecRestrictions.push_back(Cstring("4hk"));
    m_vecRestrictions.push_back(Cstring("4kl"));
    m_vecRestrictions.push_back(Cstring("4hl"));
    m_vecRestrictions.push_back(Cstring("hkl"));
    m_vecRestrictions.push_back(Cstring("W"));
    m_vecRestrictions.push_back(Cstring("X")); 
    m_vecRestrictions.push_back(Cstring("Y")); 
    m_vecRestrictions.push_back(Cstring("Z"));  
    m_vecRestrictions.push_back(Cstring("3h")); 
    m_vecRestrictions.push_back(Cstring("3l")); 
    m_vecRestrictions.push_back(Cstring("6l")); 
    m_vecRestrictions.push_back(Cstring("3hl"));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CSpaceGroupTableInfo::nGetRestrictionParity(DTREK_DWORD dwCode)
{
    std::map<DTREK_DWORD, int>::iterator         oIt;
    
    // Find restriction ID
    oIt = m_mapRestrictionParity.find(dwCode);

    if( oIt == m_mapRestrictionParity.end() )
        return -1;

    return (*oIt).second;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DTREK_DWORD CSpaceGroupTableInfo::dwGetRestrictionCode(Cstring& sName)
{    
    std::map<Cstring, DTREK_DWORD>::iterator         oIt;

    oIt = m_mapRestrictions.find(sName);

    if( oIt == m_mapRestrictions.end() )
        return DTREK_SGC_REFLN_RESTR_UNKNOWN;

    return (*oIt).second;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableInfo::bGetRestrictionInvertedLongNameByName(Cstring& sName, Cstring& sInvName)
{    
    std::map<Cstring, Cstring>::iterator         oIt;

    oIt = m_mapRestrictionInvertedLongNames.find(sName);
    if( oIt == m_mapRestrictionInvertedLongNames.end() )
    {
        sInvName = "";
        return false; // should not happen
    }

    sInvName = (*oIt).second;

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableInfo::bGetRestrictionLongNameByName(Cstring& sName, Cstring& sInvName)
{    
    std::map<Cstring, Cstring>::iterator         oIt;

    oIt = m_mapRestrictionLongNames.find(sName);
    if( oIt == m_mapRestrictionLongNames.end() )
    {
        sInvName = "";
        return false; // should not happen
    }

    sInvName = oIt->second;

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableInfo::bGetRestrictionNameByIndex(int iIndex, Cstring& sRestrName)
{
    if( iIndex < 0 || iIndex >= (int)m_vecRestrictions.size() )
        return false; // should not happen
    
    sRestrName = m_vecRestrictions[iIndex];     
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enRestr CSpaceGroupTableInfo::eGetRestrictionEnum(DTREK_DWORD dwCode)
{    
    std::map<DTREK_DWORD, enRestr>::iterator         oIt;

    oIt = m_mapRestrictionEnum.find(dwCode);
    if( oIt == m_mapRestrictionEnum.end() )
        return enRestr_Unknown;  // should not happen

    return (*oIt).second;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring& CSpaceGroupTableInfo::rsGetTableByLaueType(eLaueType eLaue)
{
    Cstring&    rs = *(&m_sTriclinic_table);
    
    switch(eLaue)
    { 
	case LAUE_1_BAR:
	    rs = *(&m_sTriclinic_table);
        break;
    case LAUE_2_OVER_M:         
	    rs = *(&m_sMonoclinic_table);
        break;
	case LAUE_MMM:      
	    rs = *(&m_sOrthorhombic_table); 
        break;
	case LAUE_4_OVER_M:
	case LAUE_4_OVER_MMM:
	    rs = *(&m_sTetragonal_table);
        break;
	case LAUE_3_BAR:
	case LAUE_3_BAR_M_1:
	case LAUE_3_BAR_1_M:
	    rs = *(&m_sTrigonal_table);
        break;
	case LAUE_6_OVER_M:
	case LAUE_6_OVER_MMM:
	    rs = *(&m_sHexagonal_table);
        break;
	case LAUE_M_3_BAR:
	case LAUE_M_3_BAR_M:
	    rs = *(&m_sCubic_table);
        break;
    }
    
    return rs;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
eLaueType CSpaceGroupTableInfo::eGetLaueClassEnumFromName(Cstring& sName)
{
    int     ii=0;
    for(ii=0; ii < (int)m_vecLaueClasses.size(); ii++)
    {
        if( m_vecLaueClasses[ii] == sName )
            break;
    }
    
    switch(ii)
    {
    case 0:
        return LAUE_1_BAR;
    case 1:
        return LAUE_2_OVER_M;
    case 2:
        return LAUE_MMM;
    case 3:
        return LAUE_4_OVER_M;
    case 4:
        return LAUE_4_OVER_MMM;
    case 5:
        return LAUE_3_BAR;
    case 6:
        return LAUE_3_BAR_M_1;
    case 7:
        return LAUE_3_BAR_1_M;
    case 8:
        return LAUE_6_OVER_M;
    case 9:
        return LAUE_6_OVER_MMM;
    case 10:
        return LAUE_M_3_BAR;
    case 11:
        return LAUE_M_3_BAR_M;
    default:
        return LAUE_UNKNOWN;
    }

    return LAUE_UNKNOWN;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring CSpaceGroupTableInfo::sGetLaueName(eLaueType eType)
{
    return m_vecLaueClasses[(int)eType];
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enCond CSpaceGroupTableInfo::eGetConditionValue(Cstring& sName)
{
    std::map<Cstring, enCond>::iterator         oIt;
    
    oIt = m_mapConditions.find(sName);

    if( oIt == m_mapConditions.end() )
        return enCond_Unknown;

    return (*oIt).second;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestrictionTableEntry::CReflectionRestrictionTableEntry()
{
    vInit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestrictionTableEntry::CReflectionRestrictionTableEntry(Cstring& sName)
{
    vInit();

    m_sName = sName;

    std::vector<Cstring>    asCodes;
    sName.nListToVector(asCodes, " ,");
    
    DTREK_DWORD       dwCode = 0L;
    
    for(int ii=0; ii < (int)asCodes.size(); ii++)
    {
        dwCode = dwGetRestrictionCode(asCodes[ii]);

        m_dwCode |= dwCode; 
        m_nParity += nGetRestrictionParity(dwCode);
        
        m_asNames.push_back(asCodes[ii]);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestrictionTableEntry::CReflectionRestrictionTableEntry(const CReflectionRestrictionTableEntry& oAnother)
{
    m_sName = oAnother.m_sName;
    m_asNames.clear();
    for(int ii=0; ii < (int)oAnother.m_asNames.size(); ii++)
    {
        m_asNames.push_back(oAnother.m_asNames[ii]);
    }
    
    m_dwCode = oAnother.m_dwCode;
    m_nParity = oAnother.m_nParity;
    
    m_enTestStatus = oAnother.m_enTestStatus;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestrictionTableEntry& 
CReflectionRestrictionTableEntry::operator=(const CReflectionRestrictionTableEntry& oAnother)
{
    m_sName = oAnother.m_sName;
    m_asNames.clear();
    for(int ii=0; ii < (int)oAnother.m_asNames.size(); ii++)
    {
        m_asNames.push_back(oAnother.m_asNames[ii]);
    }
    
    m_dwCode = oAnother.m_dwCode;
    m_nParity = oAnother.m_nParity;
    
    m_enTestStatus = oAnother.m_enTestStatus;
    
    return *this;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestrictionTableEntry::vInit()
{
    m_sName = "";
    
    m_dwCode = 0L;
    m_nParity = 0;
    
    m_enTestStatus = enUnknown;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DTREK_DWORD CReflectionRestrictionTableEntry::dwGetRestrictionCode(Cstring& sName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    return poInfo->dwGetRestrictionCode(sName);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CReflectionRestrictionTableEntry::nGetRestrictionParity(DTREK_DWORD dwCode)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    return poInfo->nGetRestrictionParity(dwCode);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionRestrictionTableEntry::bSetStatusFromRestrictionsStats
                                      (std::map<Cstring, CReflectionRestriction>& m_mapRestrictionStats)
{
    bool bFailed = false;
    bool bInconclusive  = false;

    std::map<Cstring, CReflectionRestriction>::iterator      oIt;
    enumTestStatus      eStatus;
    for(int ii=0; ii < (int)m_asNames.size(); ii++)
    {
        oIt = m_mapRestrictionStats.find(m_asNames[ii]);

        if( oIt == m_mapRestrictionStats.end() )
            return false; // should not happen

        eStatus = ((*oIt).second).eGetTestStatus();

        if( enInconclusive == eStatus )
            bInconclusive = true;
        else if( enFailed == eStatus )
            bFailed = true;
        else if( enUnknown == eStatus )
            return false; // should not happen, because all restrictions should be examined by CReflectionCondition prior to this call
    }
    
    if( bFailed )
        m_enTestStatus = enFailed;
    else if( bInconclusive )
        m_enTestStatus = enInconclusive;
    else
        m_enTestStatus = enConclusive;

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestrictionTableEntry::vGetRestrictionNames(std::vector<Cstring>& asNames)
{
    asNames.clear();
    for(int ii=0; ii < (int)m_asNames.size(); ii++)
    {
        asNames.push_back(m_asNames[ii]);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestrictionTableEntry::vSelectRestrictionsInCondition(CReflectionCondition* poCondition)
{
    for(int ii=0; ii < (int)m_asNames.size(); ii++)
        poCondition->vSelectRestriction(m_asNames[ii]);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestrictionTableEntry::vFillAbsensesArrayOldStyle(CReflectionCondition* poCondition, int nRestrictions[4])
{
    CReflectionRestriction*     poRestriction = NULL;
    
    int     ii = 0;
    for(ii=0; ii < (int)m_asNames.size(); ii++)
    {
        poRestriction = poCondition->poGetRestrictionPtrByName(m_asNames[ii]);
        poRestriction->vFillRestrictionIndexOldStyle(*&(nRestrictions[ii])); 
    }
    
    nRestrictions[ii] = enRestr_Unknown;  // to signal the end of a restrictions list for the condition
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction::CReflectionRestriction()
{
    vInit();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction::CReflectionRestriction(Cstring& sName, SPACEGROUP_CHECK_CTRL_INFO& stCtrlInfo)
{
    vInit();
    
    m_dwCode  = dwGetRestrictionCode(sName);
    m_nParity = nGetRestrictionParity(m_dwCode);

    m_sName = sName;

    m_stCtrlInfo = stCtrlInfo;
}
///////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction::CReflectionRestriction(const CReflectionRestriction& oAnother) :
CReflectionRestrictionTableEntry(oAnother)
{
    m_nHKLCompliantIntegratedReflns               = oAnother.m_nHKLCompliantIntegratedReflns;
    m_nHKLCompliantIntegratedSignificantReflns    = oAnother.m_nHKLCompliantIntegratedSignificantReflns;
    m_nHKLCompliantNonIntegratedReflns            = oAnother.m_nHKLCompliantNonIntegratedReflns;

    m_nHKLNonCompliantIntegratedReflns            = oAnother.m_nHKLNonCompliantIntegratedReflns; 
    m_nHKLNonCompliantIntegratedSignificantReflns = oAnother.m_nHKLNonCompliantIntegratedSignificantReflns; 
    m_nHKLNonCompliantNonIntegratedReflns         = oAnother.m_nHKLNonCompliantNonIntegratedReflns; 

    m_fIOverSigmaHKLNonCompliantIntegratedReflns  = oAnother.m_fIOverSigmaHKLNonCompliantIntegratedReflns; 
    m_fIOverSigmaHKLCompliantIntegratedReflns     = oAnother.m_fIOverSigmaHKLCompliantIntegratedReflns;

    m_stCtrlInfo = oAnother.m_stCtrlInfo;
    
    m_anRejects.clear();
    for(int ii=0; ii < (int)(oAnother.m_anRejects.size()); ii++)
        m_anRejects.push_back(oAnother.m_anRejects[ii]);
    
    m_bSelected = oAnother.m_bSelected;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction& CReflectionRestriction::operator=(const CReflectionRestriction& oAnother)
{
    CReflectionRestrictionTableEntry::operator=(oAnother);

    m_nHKLCompliantIntegratedReflns              = oAnother.m_nHKLCompliantIntegratedReflns;
    m_nHKLCompliantIntegratedSignificantReflns              = oAnother.m_nHKLCompliantIntegratedSignificantReflns;
    m_nHKLCompliantNonIntegratedReflns           = oAnother.m_nHKLCompliantNonIntegratedReflns;

    m_nHKLNonCompliantIntegratedReflns           = oAnother.m_nHKLNonCompliantIntegratedReflns; 
    m_nHKLNonCompliantIntegratedSignificantReflns = oAnother.m_nHKLNonCompliantIntegratedSignificantReflns; 
    m_nHKLNonCompliantNonIntegratedReflns        = oAnother.m_nHKLNonCompliantNonIntegratedReflns; 

    m_fIOverSigmaHKLNonCompliantIntegratedReflns = oAnother.m_fIOverSigmaHKLNonCompliantIntegratedReflns; 
    m_fIOverSigmaHKLCompliantIntegratedReflns    = oAnother.m_fIOverSigmaHKLCompliantIntegratedReflns;

    m_stCtrlInfo = oAnother.m_stCtrlInfo;

    m_anRejects.clear();
    for(int ii=0; ii < (int)(oAnother.m_anRejects.size()); ii++)
        m_anRejects.push_back(oAnother.m_anRejects[ii]);

    m_bSelected = oAnother.m_bSelected;

    return *this;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestriction::vInit()
{
    CReflectionRestrictionTableEntry::vInit();

    m_nHKLCompliantIntegratedReflns = 0;
    m_nHKLCompliantIntegratedSignificantReflns = 0;
    m_nHKLCompliantNonIntegratedReflns = 0;

    m_nHKLNonCompliantIntegratedReflns = 0; 
    m_nHKLNonCompliantIntegratedSignificantReflns = 0; 
    m_nHKLNonCompliantNonIntegratedReflns = 0; 

    m_fIOverSigmaHKLNonCompliantIntegratedReflns = 0.0; 
    m_fIOverSigmaHKLCompliantIntegratedReflns = 0.0;

    m_stCtrlInfo.vInit();

    m_anRejects.clear();

    m_bSelected = false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionRestriction::bTestReflection(int iIndex, int nH, int nK, int nL, double fIntensity, double fSigma)
{
    if( bIsHKLCompliant(nH, nK, nL) )
        return bTestHKLCompliantReflection(fIntensity, fSigma);
    else
        return bTestHKLNonCompliantReflection(iIndex, fIntensity, fSigma);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionRestriction::bTestHKLCompliantReflection(double fIntensity, double fSigma)
{
    if( fSigma > 0.0 )
    {
        m_nHKLCompliantIntegratedReflns++;
        m_fIOverSigmaHKLCompliantIntegratedReflns += fIntensity/fSigma;
    
        if( fIntensity/fSigma > m_stCtrlInfo.fIOverSigmaTol )
            m_nHKLCompliantIntegratedSignificantReflns++;
    
    }
    else
        m_nHKLCompliantNonIntegratedReflns++;

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionRestriction::bTestHKLNonCompliantReflection(int iIndex, double fIntensity, double fSigma)
{
    if( fSigma > 0.0 )
    {
        m_nHKLNonCompliantIntegratedReflns++;
        m_fIOverSigmaHKLNonCompliantIntegratedReflns += fIntensity/fSigma;

        if( fIntensity/fSigma > m_stCtrlInfo.fIOverSigmaTol )
        {
            m_nHKLNonCompliantIntegratedSignificantReflns++;
            
            if( m_stCtrlInfo.bCollectRejects && (int)m_anRejects.size() < m_stCtrlInfo.nRejects_MaxSaveCount )
            {
                m_anRejects.push_back(iIndex);
            }
        }    
    }
    else
        m_nHKLNonCompliantNonIntegratedReflns++;

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestriction::vTallyUpStatistics()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // First tally up I/Sigma ratios. NOTE: just because we can calculate I/Sigma here does not mean we have 
    // statistically significant number of reflections to make a conclusive determination.
    if( m_nHKLNonCompliantIntegratedReflns > 0 )
        m_fIOverSigmaHKLNonCompliantIntegratedReflns /= (double)m_nHKLNonCompliantIntegratedReflns;
    else
        m_fIOverSigmaHKLNonCompliantIntegratedReflns = 0.0;
    
    if( m_nHKLCompliantIntegratedReflns > 0 )
        m_fIOverSigmaHKLCompliantIntegratedReflns    /= (double)m_nHKLCompliantIntegratedReflns;
    else
        m_fIOverSigmaHKLCompliantIntegratedReflns = 0.0;
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // If we don't have enough integrated compliant reflections for the statistics, ther is no way
    // we can establish that non-compliant reflections are indeed absent
    if( m_nHKLCompliantIntegratedReflns < m_stCtrlInfo.nReflnCountTol_IOverSigma )
    {
        m_enTestStatus = enInconclusive;
        return;
    }
    
    // If there is an evidence that non-compliant reflections were predicted during integration,
    // we can judge whether there is a systematic absence by looking at the integration results
    // for those reflections. If, however, there is no evidence they were predicted, we must assume
    // those reflections are, indeed, absent.
    if( !(m_nHKLNonCompliantIntegratedReflns == 0 && m_nHKLNonCompliantNonIntegratedReflns == 0) )
    {
        int     nParityCoeff = m_nParity - 1;
        
        // If too few reflections have been collected/predicted we cannot judge whether there is s systematic absence or not
        if( m_nHKLNonCompliantIntegratedReflns + m_nHKLNonCompliantNonIntegratedReflns < m_stCtrlInfo.nReflnCountTol * nParityCoeff )
        {
            m_enTestStatus = enInconclusive;
            return;
        }
    }
  
    // So the remaining cases are: 
    
    // 1. Non-compliant reflections were not predicted _and_ there is enough compliant reflections
    // integrated to judge whether they are present or not. Of course, in this case we cannot calculate I/sigma for 
    // non-compliant reflections, but that is ok. We must assume they are absent (I/Sigma=0).

    // 2. Non-compliant reflections were predicted _and_ there is enough of them to judge whether they are absent or not.
    // Also, there is enough compliant reflections integrated to judge whether they are present or not. 
    // In this case we may not be able to calculate I/sigma for non-compliant reflections,
    // but that's ok - as long as there were enough of them predicted, not being able to integrate them will be 
    // considered as reflections being weak (I/sigma=0). TODO: make sure that overloads and other reasons for reflection
    // integration failure are not reported with the same error flag as weak reflections.
    
    if( m_fIOverSigmaHKLCompliantIntegratedReflns < m_stCtrlInfo.fIOverSigmaTol ||
       (m_nHKLNonCompliantIntegratedReflns >= m_stCtrlInfo.nReflnCountTol_IOverSigma && 
        m_fIOverSigmaHKLNonCompliantIntegratedReflns > m_stCtrlInfo.fIOverSigmaTol) )
    {
        m_enTestStatus = enFailed;
        return;
    }
        
    m_enTestStatus = enConclusive;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionRestriction::bIsHKLCompliant(int nH, int nK, int nL)
{
    switch(m_dwCode)
    {
    case DTREK_SGC_REFLN_RESTR_H: 
        return ((abs(nH) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_K: 
        return ((abs(nK) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_L: 
        return ((abs(nL) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_KL: 
        return ((abs(nK+nL) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_HL: 
        return ((abs(nH+nL) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_HK: 
        return ((abs(nH+nK) % 2)==0); 
    case DTREK_SGC_REFLN_RESTR_4H: 
        return ((abs(nH) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_4K: 
        return ((abs(nK) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_4L: 
        return ((abs(nL) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_4KL:
        return ((abs(nK+nL) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_4HL:
        return ((abs(nH+nL) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_4HK:
        return ((abs(nH+nK) % 4)==0);
    case DTREK_SGC_REFLN_RESTR_HKL:
        return ((abs(nH+nK+nL) % 2)==0);
    case DTREK_SGC_REFLN_RESTR_W: 
        return ((abs(-nH+nL) % 3)==0); 
    case DTREK_SGC_REFLN_RESTR_X: 
        return ((abs(nH*2+nL) % 4)==0); 
    case DTREK_SGC_REFLN_RESTR_Y: 
        return ((abs(-nH+nK+nL) % 3)==0); 
    case DTREK_SGC_REFLN_RESTR_Z: 
        return ((abs(nH-nK+nL) % 3)==0); 
    case DTREK_SGC_REFLN_RESTR_3H:
        return ((abs(nH) % 3)==0); 
    case DTREK_SGC_REFLN_RESTR_3L:
        return ((abs(nL) % 3)==0); 
    case DTREK_SGC_REFLN_RESTR_6L:
        return ((abs(nL) % 6)==0); 
    case DTREK_SGC_REFLN_RESTR_3HL:
        return ((abs(nH+nL) % 3)==0);
    default:
        return false;
    }
    
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestriction::vOutputRejects(int iTableConditionIndex,
                                            Cstring& strCondition,
                                            Creflnlist& oOriginalReflnList, 
                                            Creflnlist& oRejectsList)
{
    Cstring     sInvName("");
    vGetInvertedLongName(sInvName);

    int     iTableResrictionIndex = nGetTableRestrictionIndex();
    int     nTableCode = (nGetTableRestrictionCount() * iTableConditionIndex + iTableResrictionIndex);
    
    int         nIoverSigField    = oRejectsList.nGetFieldIndex("fIoverSig");

    int         nRestrictionField = oRejectsList.nGetFieldIndex("sRestriction");
    int         nAbsenceField     = oRejectsList.nGetFieldIndex("sAbsence");
    int         nTableCodeField   = oRejectsList.nGetFieldIndex("nTableCode");

    Crefln*     poRefln = NULL;
    int         nCurrentRejectOutputIndex = 0;

    float      fIntensity = -1.0f;
    float      fSigma = -1.0f;
    float      fIOverSigma = -1.0f;

    for(int ii=0; ii < (int)m_anRejects.size(); ii++)
    {
        poRefln = oOriginalReflnList.poGetRefln(m_anRejects[ii]);
        
        oRejectsList.nInsert(poRefln->nGetH(), 
                             poRefln->nGetK(),
                             poRefln->nGetL());
        
        nCurrentRejectOutputIndex = oRejectsList.nGetNumReflns()-1;
        
        oRejectsList[nCurrentRejectOutputIndex].vSetIntensity(fIntensity=poRefln->fGetIntensity());
        oRejectsList[nCurrentRejectOutputIndex].vSetSigmaI   (fSigma=poRefln->fGetSigmaI());
        
        fIOverSigma = fIntensity / fSigma;
        oRejectsList[nCurrentRejectOutputIndex].vSetField(nIoverSigField, fIOverSigma);
        
        oRejectsList[nCurrentRejectOutputIndex].vSetField(nRestrictionField, strCondition);
        oRejectsList[nCurrentRejectOutputIndex].vSetField(nAbsenceField,     sInvName);
        oRejectsList[nCurrentRejectOutputIndex].vSetField(nTableCodeField,   nTableCode );
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestriction::vGetInvertedLongName(Cstring& sInvName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    poInfo->bGetRestrictionInvertedLongNameByName(m_sName, sInvName);   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CReflectionRestriction::nGetTableRestrictionIndex()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    return (int)(poInfo->eGetRestrictionEnum(m_dwCode));   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CReflectionRestriction::nGetTableRestrictionCount()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    return poInfo->nGetRestrictionCount();   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionRestriction::vFillRestrictionIndexOldStyle(int& nRestrictionIndex)
{
    nRestrictionIndex = nGetTableRestrictionIndex();    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionCondition::CReflectionCondition(Cstring& sName, SPACEGROUP_CHECK_CTRL_INFO& stInfo)
{
    m_nNumberOfReflns = 0;

    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    m_eID = poInfo->eGetConditionValue(sName);
    m_strID = sName;
    
    m_stCtrlInfo = stInfo;

    m_bInconclusiveForCandidateSpaceGroup = false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB The code in bIsReflectionCondition() is mostly refactored Thad's old code in Cabsence.cpp.  
// It looks like Thad didn't finish implementing the hexagonal case.

// If we are in cubic or tetragonal we have some equivalent axes.
// Since the tables are designed to only look at absences down *one* of these
// equivalent axes, we might be ignoring important data.  In Tetragonal, we want to
// make the 0k0 and h00 restrictions match.  In Cubic, we want to make the
// h00, 0k0 and 00l restrictions match.
bool CReflectionCondition::bIsReflectionCondition(int& nH, int& nK, int& nL)
{
    bool    bCubic      = m_stCtrlInfo.eLaue == LAUE_M_3_BAR  || m_stCtrlInfo.eLaue == LAUE_M_3_BAR_M;
    bool    bTetragonal = m_stCtrlInfo.eLaue == LAUE_4_OVER_M || m_stCtrlInfo.eLaue == LAUE_4_OVER_MMM;

    switch( m_eID )
    {
    case enCond_HKL: 
        return true; // general condition
    case enCond_0KL: 
        if( nH == 0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nK==0 )
            {
	      std::swap(nH,nK);
                return true;
            }
            else if( nL==0 )
            {
                std::swap(nH,nL);
                return true;
            }
        }
        else if( bTetragonal )    
        {
            if( nK==0 ) 
            {
                std::swap(nH,nK);
                return true;
            }
        }
        return false;
    case enCond_H0L: 
        if( nK == 0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==0 )
            {
                std::swap(nH,nK);            
                return true;
            }
            else if( nL==0 )
            {
                std::swap(nL,nK);            
                return true;
            }
        }
        else if( bTetragonal )    
        {
            if( nH==0 ) 
            {
                std::swap(nH,nK);            
                return true;
            }
        }
        return false;
    case enCond_HK0: 
        if( nL == 0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==0 )
            {
                std::swap(nH,nL);            
                return true;
            }
            else if(nK==0)
            {
                std::swap(nK,nL);            
                return true;
            }
        }
        return false;
    case enCond_H00: 
        if( nK==0 && nL==0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==0 && nK==0 )
            {
                std::swap(nL,nH);
                return true;
            }
            else if( nH==0 && nL==0 ) 
            {
                std::swap(nK,nH);
                return true;
            }
        }
        else if( bTetragonal )    
        {
            if( nH==0 && nL==0 ) 
            {
                std::swap(nK,nH);
                return true;
            }
        }
        return false;
    case enCond_0K0: 
        if( nH==0 && nL==0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==0 && nK==0 )
            {
                std::swap(nK,nL);
                return true;
            }
            else if( nK==0 && nL==0 ) 
            {
                std::swap(nH,nK);            
                return true;
            }
        }
        else if( bTetragonal )    
        {
            if( nK==0 && nL==0 )
            {
                std::swap(nH,nK);
                return true;
            }
        }
        return false;
    case enCond_00L: 
        if( nH==0 && nK==0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==0 && nL==0 )
            {
                std::swap(nK,nL);
                return true;
            }
            else if(nK==0 && nL==0)    
            {
                std::swap(nH,nL);             
                return true;
            }
        }
        return false;
    case enCond_HHL: 
        if( nH==nK ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==nL )
            {
                std::swap(nK,nL);
                return true;
            }
            else if( nK==nL ) 
            {
                std::swap(nH,nL);           
                return true;
            }
        }
        return false;
    case enCond_HH0: 
        if( nH==nK && nL==0 ) // trivial case
            return true;
        else if( bCubic )
        {
            if( nH==nL && nK==0 )
            {
                std::swap(nK,nL);
                return true;
            }
            else if( nK==nL && nH==0 )
            {
                std::swap(nH,nL);               
                return true;
            }
        }
        return false;
    case enCond_HKIL: 
        return true; // general condition for hexagonal setting
    case enCond_H_H0L: 
        return (nH==-nK); 
    case enCond_HH_2_HL: 
        return (nH==nK);  
    case enCond_000L: 
        return (nH==0 && nK==0); 
    default: 
        return false; // unknown condition				
    }
    
    return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionCondition::bTestReflection(int iIndex, Crefln& oRefln)
{
    // Check if reflection is of this condition
    int     nH = oRefln.nGetH();
    int     nK = oRefln.nGetK();
    int     nL = oRefln.nGetL();
    
    //RB: just for debugging
    //if( nH == 0 || nK == 0 || nH == nK )
        //return false;
    //if( nL%4 == 0 )
        //return false;

    if( !bIsReflectionCondition(nH, nK, nL) )
        return false;
    
    m_nNumberOfReflns++;

    // Check every restriction of this condition
    std::map<Cstring, CReflectionRestriction>::iterator         oIt = m_mapRestrictions.begin();

    while( oIt != m_mapRestrictions.end() )
    {
        ((*oIt).second).bTestReflection(iIndex, nH, nK, nL,
                                       (double)oRefln.fGetIntensity(), 
                                       (double)oRefln.fGetSigmaI());

        ++oIt;
    }

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vTallyUpStatistics()
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    std::map<Cstring, CReflectionRestriction>::iterator         oIt_1 = m_mapRestrictions.begin();

    while( oIt_1 != m_mapRestrictions.end() )
    {
        ((*oIt_1).second).vTallyUpStatistics();

        ++oIt_1;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator         oIt_2 = m_mapRestrictionTableEntries.begin();

    while( oIt_2 != m_mapRestrictionTableEntries.end() )
    {
        ((*oIt_2).second).bSetStatusFromRestrictionsStats(m_mapRestrictions);

        ++oIt_2;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    vSupersedRestrictionTableEntries();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vSupersedRestrictionTableEntries()
{
    if( m_mapRestrictionTableEntries.size() < 2 )
        return; // nothing to do

    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator         oIt_1;
    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator         oIt_2;
    
    CReflectionRestrictionTableEntry*     pRest_1 = NULL;
    CReflectionRestrictionTableEntry*     pRest_2 = NULL;
    
    DTREK_DWORD     dwID_1 = 0L;
    DTREK_DWORD     dwID_2 = 0L;
    DTREK_DWORD     dwID_12 = 0L;

    int     nPar_1 = 0;
    int     nPar_2 = 0;
    
    for(oIt_1 = m_mapRestrictionTableEntries.begin(); oIt_1 != m_mapRestrictionTableEntries.end(); oIt_1++)
    {
        pRest_1 = &((*oIt_1).second);
        
        if( !pRest_1->bIsConclusive() )
            continue;
        
        dwID_1 = pRest_1->dwGetCode();
        nPar_1 = pRest_1->nGetParity();

        for(oIt_2 = m_mapRestrictionTableEntries.begin(); oIt_2 != m_mapRestrictionTableEntries.end(); oIt_2++)
        {
            pRest_2 = &((*oIt_2).second);
            dwID_2 = pRest_2->dwGetCode();
            nPar_2 = pRest_2->nGetParity();
            
            if( nPar_1 <= nPar_2 )
                continue;

            dwID_12 = dwID_1 | dwID_2;

            if( dwID_12 == dwID_1 )
                pRest_2->vSetSuperseded();
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionCondition::bAddSpaceGroupRestriction(Cstring& strRestrictionName)
{
    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator      oIt;
    oIt = m_mapRestrictionTableEntries.find(strRestrictionName);

    if( oIt != m_mapRestrictionTableEntries.end() )
        return false; // the restriction is already there!
    
    m_mapRestrictionTableEntries.insert(std::pair<Cstring, CReflectionRestrictionTableEntry>
                                       (strRestrictionName, CReflectionRestrictionTableEntry(strRestrictionName)));
    
    ////////////////////////////////////////////////////////////////
    std::vector<Cstring>        asNames;
    oIt = m_mapRestrictionTableEntries.find(strRestrictionName);
    ((*oIt).second).vGetRestrictionNames(asNames);
    ///////////////////////////////////////////////////////////////
    
    for(int ii=0; ii < (int)asNames.size(); ii++)
        bAddRestriction(asNames[ii]);

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CReflectionCondition::bAddRestriction(Cstring& strRestrictionName)
{
    std::map<Cstring, CReflectionRestriction>::iterator      oIt;
    oIt = m_mapRestrictions.find(strRestrictionName);
    
    if( oIt != m_mapRestrictions.end() )
        return false; // the restriction is already there!
    
    m_mapRestrictions.insert(std::pair<Cstring, CReflectionRestriction>
                            (strRestrictionName, CReflectionRestriction(strRestrictionName, m_stCtrlInfo)));

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vOutputRejects(Creflnlist& oOriginalReflnList, Creflnlist& oRejectsList)
{
    int     iConditionTableIndex = (int)m_eID;

    std::map<Cstring, CReflectionRestriction>::iterator         oIt = m_mapRestrictions.begin();

    while( oIt != m_mapRestrictions.end() )
    {
        ((*oIt).second).vOutputRejects(iConditionTableIndex, m_strID, oOriginalReflnList, oRejectsList);

        ++oIt;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vMakeFullRestrictionList()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    std::map<Cstring, DTREK_DWORD>    mapRestrictionNamesIDs = poInfo->rmapGetRestrictionNamesIDs();
    std::map<Cstring, DTREK_DWORD>::iterator         oIt;
    
    Cstring     sTemp("");
    Cstring&    rs = sTemp;

    for(oIt=mapRestrictionNamesIDs.begin(); oIt != mapRestrictionNamesIDs.end(); oIt++)
    {
        rs = *&(oIt->first);
        bAddRestriction(rs);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enumTestStatus CReflectionCondition::eGetRestrictionTableEntryStatus(Cstring& strRestrictionName)
{
    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator       oIt;

    oIt = m_mapRestrictionTableEntries.find(strRestrictionName);
    
    if( oIt == m_mapRestrictionTableEntries.end() )
        return enUnknown; // should not happen

    return ((*oIt).second).eGetTestStatus();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction*  CReflectionCondition::poGetRestrictionPtrByName(Cstring& sName)
{
    std::map<Cstring, CReflectionRestriction>::iterator       oIt;
    oIt = m_mapRestrictions.find(sName);

    if( oIt == m_mapRestrictions.end() )
        return NULL;

    return &(oIt->second);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestrictionTableEntry*  CReflectionCondition::poGetRestrictionTableEntryPtrByName(Cstring& sName)
{
    std::map<Cstring, CReflectionRestrictionTableEntry>::iterator       oIt;
    oIt = m_mapRestrictionTableEntries.find(sName);

    if( oIt == m_mapRestrictionTableEntries.end() )
        return NULL;

    return &(oIt->second);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vSelectRestriction(Cstring& sName)
{
    CReflectionRestriction*     poRestr = poGetRestrictionPtrByName(sName); 
    
    if( !poRestr )
        return;    // Shouldn't happen. Just a safety check

    poRestr->vSelect(true);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CReflectionCondition::vUnselectAllRestrictions()
{
    std::map<Cstring, CReflectionRestriction>::iterator     oIt;
    for(oIt=m_mapRestrictions.begin(); oIt != m_mapRestrictions.end(); oIt++)
    {
        (oIt->second).vSelect(false);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CReflectionCondition::nGetNumberOfSelectedRestrictions()
{
    int nCount = 0;

    std::map<Cstring, CReflectionRestriction>::iterator     oIt;
    for(oIt=m_mapRestrictions.begin(); oIt != m_mapRestrictions.end(); oIt++)
    {
        if( (oIt->second).bIsSelected() )
            nCount++;
    }
    
    return nCount;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionRestriction* CReflectionCondition::poGetFirstSelectedRestrictionPtr()
{
    std::map<Cstring, CReflectionRestriction>::iterator     oIt;
    for(oIt=m_mapRestrictions.begin(); oIt != m_mapRestrictions.end(); oIt++)
    {
        if( (oIt->second).bIsSelected() )
            return &(oIt->second);
    }
    
    return NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CSpaceGroupCheckSolution::CSpaceGroupCheckSolution(int iIndex, 
                                                   int nNumber, 
                                                   Cstring& sName,
                                                   eLaueType eLaue,
                                                   bool bCentric, 
                                                   int nFreq,
                                                   bool bNonStandard)
{
    m_nIndexInLaueClassTable = iIndex;
    m_nIntTablesNumber = nNumber;
    
    m_eLaue = eLaue;
    m_bCentric = bCentric;
    m_bChiral = bIsChiralCompatibleSpacegroup(nNumber);
    
    m_sName = sName;
    
    ///////////////////////////////////////////////////
    // Figure out presentation
    // The difference between m_sName and m_sPresentation is
    // that m_sName comes directly from Thad's spacegroup check table, so m_sName 
    // for some monoclinic cells has cell choice information in parentheses, e.g. (C2), (C3). 
    // m_sPresentation is just a spacegroup name without cell choice information.
    
    m_sPresentation = sName.before('(');

    m_nFrequency = nFreq;
    m_fProbability = 0.0;

    m_bNonStandard = bNonStandard;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CSpaceGroupCheckSolution::CSpaceGroupCheckSolution(const CSpaceGroupCheckSolution& oAnother)
{
    m_nIndexInLaueClassTable = oAnother.m_nIndexInLaueClassTable;
    m_nIntTablesNumber       = oAnother.m_nIntTablesNumber;
    
    m_bCentric               = oAnother.m_bCentric;
    m_eLaue                  = oAnother.m_eLaue;
    m_bChiral                = oAnother.m_bChiral;
    
    m_sName                  = oAnother.m_sName;
    m_sPresentation          = oAnother.m_sPresentation;

    m_nFrequency             = oAnother.m_nFrequency;
    m_fProbability           = oAnother.m_fProbability;

    m_bNonStandard           = oAnother.m_bNonStandard;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If a space group contains mirror planes or inversion center, it cannot be chiral
bool CSpaceGroupCheckSolution::bIsChiralCompatibleSpacegroup(int nSpaceGroup_IntTables)
{
    char*   pcName = &(g_cpSpaceGroupNames[nSpaceGroup_IntTables - 1][1]);
    
    return (strspn(pcName, "0123456789") == strlen(pcName));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring CSpaceGroupCheckSolution::sGetIntTablesName()
{
    return Cstring(g_cpSpaceGroupNames[m_nIntTablesNumber - 1]);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableEntry::bTest(std::vector<CReflectionCondition*>& vecConditionPtrs)
{
    m_strAbsenceTest = "";
    
    std::map<Cstring, Cstring>::iterator     oIt;

    bool    bAtLeastOneInconclusiveTest = false;
    bool    bAtLeastOneSupersededTest  = false;
    bool    bAtLeastOneFailedTest  = false;
    
    enumTestStatus  eStatus = enUnknown;

    for(int ii=0; ii < (int)vecConditionPtrs.size(); ii++)
    {
        oIt = m_mapConditionsRestrictions.find(vecConditionPtrs[ii]->m_strID);

        if( oIt == m_mapConditionsRestrictions.end() )
        {
            m_strAbsenceTest += '-';
            continue;
        }
        
        eStatus = vecConditionPtrs[ii]->eGetRestrictionTableEntryStatus(*&(oIt->second));

        if( eStatus == enUnknown )
            return false;   // should not happen!
        else if( eStatus == enInconclusive )
        {
            m_strAbsenceTest += '?';
            bAtLeastOneInconclusiveTest = true;
        }
        else if( eStatus == enSuperseded )
        {
            m_strAbsenceTest += 'S';
            bAtLeastOneSupersededTest = true;
        }
        else if( eStatus == enFailed )
        {
            m_strAbsenceTest += 'F';
            bAtLeastOneFailedTest = true;
        }
        else
        {
            m_strAbsenceTest += 'A';
            m_nConclusivePassedTests++;
        }
    }
        
    if( bAtLeastOneFailedTest )
        m_eTestStatus = enFailed;
    else if( bAtLeastOneSupersededTest )
        m_eTestStatus = enSuperseded;
    else if( bAtLeastOneInconclusiveTest )
        m_eTestStatus = enInconclusive;
    else 
        m_eTestStatus = enConclusive;

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CSpaceGroupTableEntry::vMarkInconclusiveCondtions(std::vector<CReflectionCondition*>& vecConditionPtrs)
{
#ifdef _DEBUG
    if( m_strAbsenceTest.length() != vecConditionPtrs.size() )
    {
        printf("\n\n\n!!!!!!!!! ERROR !!!!!!!!!\n\n\n"); // shouldn't happen
        
        return;
    }
#endif

    for(int ii=0; ii < (int)m_strAbsenceTest.GetLength(); ii++)
    {
        if( m_strAbsenceTest.GetAt(ii) == '?' )
            vecConditionPtrs[ii]->vSetInconclusiveForCandidateSpaceGroup();
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableEntry::bGenerateSpaceGroupSolutions(Cstring& sLaueClassAndCentricityInfo, 
                                                         Cstring& sSpaceGroupsInfo)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine Laue class
    std::vector<Cstring>    vecLaueClassAndCentricityInfo;
    sLaueClassAndCentricityInfo.nListToVector(vecLaueClassAndCentricityInfo, ".");

    if( vecLaueClassAndCentricityInfo.size() != 2 )
        return false; // should not happen, because we must always have "Laue class" and "Centricity",
                      // that is exactly two pieces of information
    
    eLaueType    eLaue = eGetLaueClassEnumFromName(vecLaueClassAndCentricityInfo[0]);
    if( eLaue == LAUE_UNKNOWN )
        return false; // should not happen, because Laue class must be known
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine Centricity
    SPACEGROUP_CHECK_CTRL_INFO::eCentricityType     eCent = SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Unknown;

    if( vecLaueClassAndCentricityInfo[1].contains("AC") )
        eCent = SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric;
    else if( vecLaueClassAndCentricityInfo[1].contains("C") )
        eCent = SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric;
    else
        return false; // should not happen
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine space group numbers, names, etc
    char        pcBuf[80];
    strcpy(pcBuf, sSpaceGroupsInfo.string());
    char*       pc=&pcBuf[0];
    
    char*       pc2 = NULL;
    char*       pc3 = NULL;

    Cstring     sSpaceGroupName("");
    int         nSpaceGroup_IntTables = -1;
    
    // Load the space group names, and the space group numbers.
    bool    bNonStandard = false;
    while( pc[0] )
    {
        if( pc2 = strchr(pc, ',') ) 
        { 
            *pc2=0;
            
            std::swap(pc2,pc);
            
            pc++; 
        } 
        else 
        { 
            pc2=pc; 
            pc+=strlen(pc); 
        }
        
        if( strcmp(pc2,"*") ) // if pc2 is not pointing to '*'
        {
            // Look for a '.'.  This indicates that a non-standard space group will need to be changed. (if we are reindexing).
            if( pc3 = strchr(pc2,'.') )
            {
                ///////////////////////////////////////////////
                //strcpy(cpSpaceNames[nSpaceGroupCt],pc3+1);
                sSpaceGroupName = pc3 + 1;
                ///////////////////////////////////////////////

                *pc3=0;
                
                //if( 1 != sscanf(pc2, "%d", &nSpaceGroups[nSpaceGroupCt]) )
                    //ERROR;
                
                nSpaceGroup_IntTables = atoi(pc2);

                bNonStandard = true;
            }
            else
            {
                //if( 1 != sscanf(pc2,"%d",&nSpaceGroups[nSpaceGroupCt]) )
                    //ERROR;
                nSpaceGroup_IntTables = atoi(pc2);

                //strcpy(cpSpaceNames[nSpaceGroupCt], g_cpSpaceGroupNames[nSpaceGroups[nSpaceGroupCt]-1]);
                sSpaceGroupName = g_cpSpaceGroupNames[nSpaceGroup_IntTables - 1];
                
                bNonStandard = false;
            }
            
            m_vecSpaceGroupSolutions.push_back(CSpaceGroupCheckSolution(m_nIndexInLaueClassTable, 
                                                                        nSpaceGroup_IntTables,
                                                                        sSpaceGroupName,
                                                                        eLaue,
                                                                        eCent == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric ?
                                                                          true : false,
                                                                        CSpaceGroupTableInfo::m_nSpaceGroupFreq[nSpaceGroup_IntTables-1],
                                                                        bNonStandard));
        }
    }
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
eLaueType CSpaceGroupTableEntry::eGetLaueClassEnumFromName(Cstring& sName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    return poInfo->eGetLaueClassEnumFromName(sName);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The idea of this function is that if a space group table entry is selected, then all restrictions 
// referenced in that entry need to be selected for output. The problem is that restriction objects are 
// not part of a space group table entry object. Instead, they are part of condition objects. So we
// need to find the right condition objects and select the right restriction objects in them.
// We have to start with "restriction table entry" objects, and using the names of restrictions, referenced in them,
// make our way to the restriction objects and select them.
void CSpaceGroupTableEntry::vSelectRestrictionsInConditions(std::vector<CReflectionCondition*>& vecConditionPtrs)
{
    Cstring     sTemp("");
    
    std::map<Cstring, Cstring>::iterator        oIt;
    
    CReflectionRestrictionTableEntry*    poRestrEntry = NULL;

    for(int ii=0; ii < (int)vecConditionPtrs.size(); ii++)
    {
        sTemp = vecConditionPtrs[ii]->sGetID();
        
        if( (oIt = m_mapConditionsRestrictions.find(sTemp)) != m_mapConditionsRestrictions.end() )
        {
            sTemp = oIt->second; // The name of a composite restriction for this condition    
            
            poRestrEntry = vecConditionPtrs[ii]->poGetRestrictionTableEntryPtrByName(sTemp);
            
            poRestrEntry->vSelectRestrictionsInCondition(vecConditionPtrs[ii]);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enumTestStatus CSpaceGroupTableEntry::eGetConditionStatus(CReflectionCondition* poCond)
{
    std::map<Cstring, Cstring>::iterator        oIt = m_mapConditionsRestrictions.find(poCond->sGetID());
    
    if( oIt == m_mapConditionsRestrictions.end() )
        return enUnknown;

    return poCond->eGetRestrictionTableEntryStatus(oIt->second);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring CSpaceGroupTableEntry::sGetConditionRestriction(CReflectionCondition* poCond)
{
    std::map<Cstring, Cstring>::iterator        oIt = m_mapConditionsRestrictions.find(poCond->sGetID());
    
    if( oIt == m_mapConditionsRestrictions.end() )
        return Cstring("");

    return (oIt->second);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CSpaceGroupTableEntry::vFillAbsensesArrayOldStyle(std::vector<CReflectionCondition*>& vecConditionPtrs, 
                                                       int nAbsences[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4])
{
    std::map<Cstring, Cstring>::iterator        oIt;
    Cstring     sRestrictionTableEntryname("");
    CReflectionRestrictionTableEntry*           poRestrTableEntry = NULL;

    for(int ii=0; ii < (int)vecConditionPtrs.size(); ii++)
    {
        oIt = m_mapConditionsRestrictions.find(vecConditionPtrs[ii]->sGetID());
        
        if( oIt == m_mapConditionsRestrictions.end() )
        {
            nAbsences[ii][0] = enRestr_Unknown;  // to signal the end of a restrictions list for the condition
            continue;
        }
        
        sRestrictionTableEntryname = oIt->second;
        
        poRestrTableEntry = vecConditionPtrs[ii]->poGetRestrictionTableEntryPtrByName(sRestrictionTableEntryname);

        poRestrTableEntry->vFillAbsensesArrayOldStyle(vecConditionPtrs[ii], nAbsences[ii]);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return number of added solutions
int CSpaceGroupTableEntry::nAddSpaceGroupSolutionsToVector(std::vector<CSpaceGroupCheckSolution>& vecSpaceGroupSolutions,
                                                           eLaueType eLaue,
                                                           bool bCentric,
                                                           bool bSelectChiralOnly)
{
    int     nRet = 0;
    
    for(int ii=0; ii < (int)m_vecSpaceGroupSolutions.size(); ii++)
    {
        if( eLaue != LAUE_UNKNOWN && eLaue != m_vecSpaceGroupSolutions[ii].eGetLaueType() )
            continue;

        if( m_vecSpaceGroupSolutions[ii].bIsCentric() && !bCentric )
            continue;

        if( !m_vecSpaceGroupSolutions[ii].bIsCentric() && bCentric )
            continue;

        if( !m_vecSpaceGroupSolutions[ii].bIsChiral() && bSelectChiralOnly )
            continue;

        vecSpaceGroupSolutions.push_back(m_vecSpaceGroupSolutions[ii]);
        nRet++;
    }

    return nRet;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CSpaceGroupTableEntry::bIsSpaceGroupInStandardSetting(int nSG)
{
    for(int ii=0; ii < (int)m_vecSpaceGroupSolutions.size(); ii++)
    {
        if( nSG == m_vecSpaceGroupSolutions[ii].nGetIntTablesNumber() &&
            m_vecSpaceGroupSolutions[ii].bIsStandard() )
            return true;
    }
    
    return false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cspacegroup_check::Cspacegroup_check(SPACEGROUP_CHECK_CTRL_INFO& stInfo)
{
    m_stCtrlInfo = stInfo;
    
    m_nSelectedCandidateSpaceGroupIndex = -1;

    vLoadTableByLaueType();

    vMakeFullRestrictionListForConditionHKL();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cspacegroup_check::~Cspacegroup_check()
{
    int     ii = 0;
    for(ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
        delete m_vecReflnConditionPtrs[ii];

    m_vecReflnConditionPtrs.clear();

    for(ii=0; ii < (int)m_vecSpaceGroupTableEntryPtrs.size(); ii++)
        delete m_vecSpaceGroupTableEntryPtrs[ii];
    
    m_vecSpaceGroupTableEntryPtrs.clear();
    
    DoSpacegroupTableInfo(DTREK_DSGTI_DELETE);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cspacegroup_check::bTestReflnList(Creflnlist& oReflnList)
{
    Crefln*     poRefln = NULL;
    int         ii = 0;
    for(int nRef=0; nRef < oReflnList.nGetNumReflns(); nRef++)
    {
		poRefln = oReflnList.poGetRefln(nRef);
        
        for(ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
        {
            m_vecReflnConditionPtrs[ii]->bTestReflection(nRef, *poRefln);
        }
        
        vDoTexsanStats(*poRefln);
    }
    /////////////////////////////////////////////////////////////////////
    for(ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        m_vecReflnConditionPtrs[ii]->vTallyUpStatistics();
    }
    ////////////////////////////////////////////////////////////////////
    

    m_stTexsanStats.vTally();
    
    /////////////////////////////////////////////////////////////////////
    if( m_stCtrlInfo.bCollectRejects )
        vOutputRejects(oReflnList);

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vLoadTableByLaueType()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    Cstring&    roTable = poInfo->rsGetTableByLaueType(m_stCtrlInfo.eLaue);
    
    vLoadTable(roTable);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vLoadTable(Cstring& sInput)
{
    // Define white space
    char* pcWhiteSpace = "\t\n\r ";
    
    std::vector<Cstring>        asTokens;
    sInput.nListToVector(asTokens, pcWhiteSpace);
    
    int     nTotalNumberOfTokens = (int)asTokens.size();
    
    CSpaceGroupTableEntry*      pSpaceGroupTableEntry = NULL;

    // First read in the conditions
    int     ii = 1; // The very first token (ii=0) must be a line offset, so we are just skipping it here
    while( asTokens[ii] != "|" )
    {
        vAddCondition(asTokens[ii]);
        ii++;
    }
    
    std::vector<Cstring>    vecCentricAcentricTitles;
    ii++; // skip "|"
    while( asTokens[ii] != "STOP" )
    {
        vecCentricAcentricTitles.push_back(asTokens[ii]);
        ii++;
    }

    int     nNumberOfConditions = (int)m_vecReflnConditionPtrs.size();
    
    ii += 2; // skip STOP and the first token after STOP
    int     nSpaceGroupEntryIndex = -1;
    
    Cstring     sTemp("");
    while( ii < nTotalNumberOfTokens )
    {
        pSpaceGroupTableEntry = new CSpaceGroupTableEntry(++nSpaceGroupEntryIndex);
        m_vecSpaceGroupTableEntryPtrs.push_back(pSpaceGroupTableEntry);

        for(int jj=0; jj < nNumberOfConditions; jj++, ii++)
        {
            if( asTokens[ii] != "-" )
            {
                bAddRestrictionToCondition(jj, asTokens[ii]);
                
                sTemp = m_vecReflnConditionPtrs[jj]->sGetID();
                std::pair<Cstring, Cstring>     pair(sTemp, asTokens[ii]);
                pSpaceGroupTableEntry->vAddConditionRestriction(pair);
            }
        }
        
        for(int kk=0; kk < (int)vecCentricAcentricTitles.size(); kk++, ii++)
        {
            pSpaceGroupTableEntry->bGenerateSpaceGroupSolutions(vecCentricAcentricTitles[kk], asTokens[ii]);
        }
        
        ii++; // now we have moved over to the next spacegroup table entry first condition column
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vAddCondition(Cstring& sCondition)
{
    CReflectionCondition*       pCond = new CReflectionCondition(sCondition, m_stCtrlInfo);
    m_vecReflnConditionPtrs.push_back(pCond);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CReflectionCondition* Cspacegroup_check::poGetConditionByID(enCond eID)
{
    for(int ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        if( eID == m_vecReflnConditionPtrs[ii]->eGetID() )
            return m_vecReflnConditionPtrs[ii];
    }
    
    return NULL;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vMakeFullRestrictionListForConditionHKL()
{
    CReflectionCondition*   poConditionHKL = poGetConditionByID(enCond_HKL);
    
    Cstring     sTemp("hkl");
    
    if( !poConditionHKL )
    {
        vAddCondition(sTemp);
        poConditionHKL = poGetConditionByID(enCond_HKL);
    }    
    
    poConditionHKL->vMakeFullRestrictionList();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cspacegroup_check::bAddRestrictionToCondition(int iConditionIndex, 
                                                   Cstring& strRestrictionName)
{
    m_vecReflnConditionPtrs[iConditionIndex]->bAddSpaceGroupRestriction(strRestrictionName);
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vOutputRejects(Creflnlist& oOriginalReflnList)
{
    Creflnlist      oRejectsList;

    int     nRestriction    = oRejectsList.nExpandGetField("sRestriction");  // RB: this is really called "condition"
    int     nAbsence        = oRejectsList.nExpandGetField("sAbsence");      // RB: this is really called "restriction"
    int     nTableCode      = oRejectsList.nExpandGetField("nTableCode");
    int     nIoverSig       = oRejectsList.nExpandGetField("fIoverSig");
    
    for(int ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        m_vecReflnConditionPtrs[ii]->vOutputRejects(oOriginalReflnList, oRejectsList);
    }
    
    oRejectsList.vSort(eReflnField_int_type, nTableCode, NULL);

    printf("Outputing 'dtcell_rejects.ref'...\n");
    
    oRejectsList.nWrite("dtcell_rejects.ref");
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vDoTexsanStats(Crefln& oRefln)
{
    int nH = oRefln.nGetH();
    int nK = oRefln.nGetK();
    int nL = oRefln.nGetL();
    
    double  fIntensity = (double)oRefln.fGetIntensity();
    double  fSigma     = (double)oRefln.fGetSigmaI();
    
    double  fIoverSigma = 0.0;
    
    if( fSigma > 0.0 )
        fIoverSigma = fIntensity / fSigma;
    
    bool bAbsent = fIoverSigma < m_stCtrlInfo.fIOverSigmaTol;
    
    int nHKLParity = 0;
    int nHKParity  = 0;
    int nHLParity  = 0;
    int nKLParity  = 0;

    // Discover parity of the HKL indices.
    if (abs(nL) % 2) 
    {
        nHKLParity+=1;
        nHLParity+=1;
        nKLParity+=1;
    }

    if (abs(nK) % 2)
    {
        nHKLParity+=2;
        nHKParity+=1;
        nKLParity+=2;
    }

    if (abs(nH) % 2)
    {
        nHKLParity+=4;
        nHKParity+=2;
        nHLParity+=2;
    }

    m_stTexsanStats.vAdd(nHKLParity          + TEX_EEE, bAbsent, fIoverSigma);
    
    int     nAdd = -1;

    if ( nL == 0 /*abRestrictions[RES_HK0]*/ )
        m_stTexsanStats.vAdd(nHKParity       + TEX_HK0_EE, bAbsent, fIoverSigma);
    
    if ( nK == 0 /*abRestrictions[RES_H0L]*/ )
        m_stTexsanStats.vAdd(nHLParity       + TEX_H0L_EE, bAbsent, fIoverSigma);
    
    if ( nH == 0 /*abRestrictions[RES_0KL]*/)
        m_stTexsanStats.vAdd(nKLParity       + TEX_0KL_EE, bAbsent, fIoverSigma);
    
    if ( nH == nK /*abRestrictions[RES_HHL]*/)
        m_stTexsanStats.vAdd(nHLParity       + TEX_HHL_EE, bAbsent, fIoverSigma);
    
    if ( nH == -nK )
        m_stTexsanStats.vAdd(nHLParity       + TEX_H_HL_EE, bAbsent, fIoverSigma);
    
    if ( nH == nK && nK == nL )
        m_stTexsanStats.vAdd((abs(nH) % 2)   + TEX_HHH_E, bAbsent, fIoverSigma);
    
    if ( nH == nK && nL == 0 )
        m_stTexsanStats.vAdd((abs(nH) % 2)   + TEX_HH0_E, bAbsent, fIoverSigma);
    
    if ( nK == 0 && nL == 0 /*abRestrictions[RES_H00]*/)
        m_stTexsanStats.vAdd((abs(nH) % 2)   + TEX_H00_E, bAbsent, fIoverSigma);
    
    if ( nH == 0 && nL == 0 /*abRestrictions[RES_0K0]*/)
        m_stTexsanStats.vAdd((abs(nK) % 2)   + TEX_0K0_E, bAbsent, fIoverSigma);
    
    if ( nH == 0 && nK == 0 /*abRestrictions[RES_00L]*/)
        m_stTexsanStats.vAdd((abs(nL) % 2)   + TEX_00L_E, bAbsent, fIoverSigma);
    
    if ( nH == 0 /*abRestrictions[RES_0KL]*/) 
    {
        nAdd = abs(nK + nL)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd            + TEX_0KL_A /*((abAbsences[ABS_4KL])?(0):(1))*/, bAbsent, fIoverSigma);
    }    
    
    if ( nK == 0 /*abRestrictions[RES_H0L]*/) 
    {
        nAdd = abs(nH + nL)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_H0L_A /*((abAbsences[ABS_4HL])?(0):(1))*/, bAbsent, fIoverSigma);
    }

    if ( nL == 0 /*abRestrictions[RES_HK0]*/) 
    {
        nAdd = abs(nH + nK)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_HK0_A /*((abAbsences[ABS_4HK])?(0):(1))*/, bAbsent, fIoverSigma);
    }

    if ( nH == 0 && nL == 0 /*abRestrictions[RES_0K0]*/) 
    {
        nAdd = abs(nK)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_0K0_A /*((abAbsences[ABS_4K])?(0):(1))*/  , bAbsent, fIoverSigma);
    }

    if ( nK == 0 && nL == 0 /*abRestrictions[RES_H00]*/) 
    {
        nAdd = abs(nH)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_H00_A /*((abAbsences[ABS_4H])?(0):(1))*/  , bAbsent, fIoverSigma);
    }
    
    if ( nH == 0 && nK == 0 /*abRestrictions[RES_00L]*/) 
    {
        nAdd = abs(nL)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_00L_A /*((abAbsences[ABS_4H])?(0):(1))*/  , bAbsent, fIoverSigma);
    }    
    
    if ( nH == nK /*abRestrictions[RES_HHL]*/) 
    {
        nAdd = abs(2*nH + nL)%4 == 0 ? 0 : 1;
        m_stTexsanStats.vAdd(nAdd           + TEX_HHL_A/*((abAbsences[ABS_T])?(0):(1))*/    , bAbsent, fIoverSigma);
    }

    if ( nH == -nK /*abRestrictions[RES_H_H0L]*/) 
    {
        if (((abs(nH+nL) % 3) == 0) && (abs(nL) % 2))
            m_stTexsanStats.vAdd(TEX_H_H0L_A1, bAbsent, fIoverSigma);
        if ((abs(nH+nL) % 3) == 0)
            m_stTexsanStats.vAdd(TEX_H_H0L_A2, bAbsent, fIoverSigma);
        if ((abs(nH+nL) % 3) != 0)
            m_stTexsanStats.vAdd(TEX_H_H0L_A3, bAbsent, fIoverSigma);
        if (((abs(-nH+nL) % 3) == 0) && ((abs(nL) % 2)!=0))
            m_stTexsanStats.vAdd(TEX_H_H0L_B1, bAbsent, fIoverSigma);
        if ((abs(-nH+nL) % 3) == 0)
            m_stTexsanStats.vAdd(TEX_H_H0L_B2, bAbsent, fIoverSigma);
        if ((abs(-nH+nL) % 3) != 0)
            m_stTexsanStats.vAdd(TEX_H_H0L_B3, bAbsent, fIoverSigma);
    }
    
    if( nH == 0 && nK == 0 /*abRestrictions[RES_000L]*/)
    {
        if ((abs(nL) % 3) == 0)
            m_stTexsanStats.vAdd(TEX_000L_A, bAbsent, fIoverSigma);
        if ((abs(nL) % 3) != 0)
            m_stTexsanStats.vAdd(TEX_000L_NA, bAbsent, fIoverSigma);
        if ((abs(nL) % 6) == 0)
            m_stTexsanStats.vAdd(TEX_000L_B, bAbsent, fIoverSigma);
        if ((abs(nL) % 6) != 0)
            m_stTexsanStats.vAdd(TEX_000L_NB, bAbsent, fIoverSigma);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vTestSpaceGroupTableEntries()
{
    for(int ii=0; ii < (int)m_vecSpaceGroupTableEntryPtrs.size(); ii++)
    {
        m_vecSpaceGroupTableEntryPtrs[ii]->bTest(m_vecReflnConditionPtrs);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//mrpvoid Cspacegroup_check::vSelectCandidateSpaceGroups()
void Cspacegroup_check::vSelectCandidateSpaceGroups(const int nFavoredSpaceGroupNumber) //mrp
{    
    int     ii = 0;
    int     nMaxNumberOfConclusiveTests = 0;
    // First find out if we have at least one conclusive entry except for the very first one,
    // which does not have any restrictions i.e. tests to perform
    for(ii=1; ii < (int)m_vecSpaceGroupTableEntryPtrs.size(); ii++)
    {
        if( m_vecSpaceGroupTableEntryPtrs[ii]->bIsConclusive() )
        {
            if( m_vecSpaceGroupTableEntryPtrs[ii]->nGetNumberOfConclusivePassedTests() > nMaxNumberOfConclusiveTests )
                nMaxNumberOfConclusiveTests = m_vecSpaceGroupTableEntryPtrs[ii]->nGetNumberOfConclusivePassedTests();
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<CSpaceGroupTableEntry*>     vecCandidateSpaceGroupTableEntryPtrs;
    
    ii = nMaxNumberOfConclusiveTests > 0 ? 1 : 0; // if there is at least one table entry that conclusively passed the test, 
                                                  // then we do not need the very first table entry (which does not have any
                                                  // restrictions i.e. tests to perform). 

    for(; ii < (int)m_vecSpaceGroupTableEntryPtrs.size(); ii++)
    {
        if( ii == 0 ||  // if we started from index 0, it means we don't have any conclusive tests besides entry with index 0
            m_vecSpaceGroupTableEntryPtrs[ii]->bIsInconclusive() ||
            m_vecSpaceGroupTableEntryPtrs[ii]->bIsConclusive() && m_vecSpaceGroupTableEntryPtrs[ii]->nGetNumberOfConclusivePassedTests() 
                                                                  == nMaxNumberOfConclusiveTests
          )
        {
            m_vecSpaceGroupTableEntryPtrs[ii]->vSetSelected();
            vecCandidateSpaceGroupTableEntryPtrs.push_back(m_vecSpaceGroupTableEntryPtrs[ii]);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if( vecCandidateSpaceGroupTableEntryPtrs.size() > 1 )
        vDeselectCandidateSpaceGroupTableEntriesByComparison(vecCandidateSpaceGroupTableEntryPtrs);
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Now generate spacegroup solutions based on the selected spacegroup table entries as well as on the user input
    // for Laue class, centricity, chirality. Additionally, if a table entry does not result in any solutions - deselect it.
    int     nSolutionCountPerTableEntry = 0;
    for(ii=0; ii < (int)vecCandidateSpaceGroupTableEntryPtrs.size(); ii++)
    {
        nSolutionCountPerTableEntry = 0;

        if( !vecCandidateSpaceGroupTableEntryPtrs[ii]->bIsSelected() )
            continue;

        nSolutionCountPerTableEntry += nGetSpaceGroupSolutionsFromSpaceGroupTableEntry(vecCandidateSpaceGroupTableEntryPtrs[ii]->nGetIndexInLaueClassTable(), 
                                                                                       SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric);

        nSolutionCountPerTableEntry += nGetSpaceGroupSolutionsFromSpaceGroupTableEntry(vecCandidateSpaceGroupTableEntryPtrs[ii]->nGetIndexInLaueClassTable(),
                                                                                       SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric);
        if( 0 == nSolutionCountPerTableEntry )
            vecCandidateSpaceGroupTableEntryPtrs[ii]->vSetSelected(false);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vMarkInconclusiveConditionsForCandidateSpaceGroupTableEntries(vecCandidateSpaceGroupTableEntryPtrs);
    
    vFindProbabilitiesForCandidateSpaceGroups();

    int     nBestCentricSpaceGroupIndex  = nFindBestCandidateSpaceGroup(SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric);
    int     nBestAcentricSpaceGroupIndex = nFindBestCandidateSpaceGroup(SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric);
    
    m_nSelectedCandidateSpaceGroupIndex = -1;
    
    if( m_stCtrlInfo.eCentricity == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric )
    {
        // Give preference to centric group
        if( nBestCentricSpaceGroupIndex > -1 )
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestCentricSpaceGroupIndex;
        }
        else if ( nBestAcentricSpaceGroupIndex > -1 )
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestAcentricSpaceGroupIndex; 
        }
    }
    else if( m_stCtrlInfo.eCentricity == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric )
    {
        // Give preference to acentric group
        if( nBestAcentricSpaceGroupIndex > -1 )
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestAcentricSpaceGroupIndex;
        }
        else if ( nBestCentricSpaceGroupIndex > -1 )
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestCentricSpaceGroupIndex; 
        }
    }
    else if( m_stCtrlInfo.eCentricity == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Unknown )
    {
        if( nBestCentricSpaceGroupIndex > -1 && nBestAcentricSpaceGroupIndex == -1 ) 
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestCentricSpaceGroupIndex; 
        }
        else if( nBestAcentricSpaceGroupIndex > -1 && nBestCentricSpaceGroupIndex == -1  )    
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestAcentricSpaceGroupIndex; 
        }
        else if( nBestAcentricSpaceGroupIndex > -1 && nBestCentricSpaceGroupIndex > -1 )
        {
            m_nSelectedCandidateSpaceGroupIndex = nBestAcentricSpaceGroupIndex; 

            if( m_vecCandidateSpaceGroups[nBestCentricSpaceGroupIndex].fGetProbability() > 
                m_vecCandidateSpaceGroups[nBestAcentricSpaceGroupIndex].fGetProbability() )
            {
                m_nSelectedCandidateSpaceGroupIndex = nBestCentricSpaceGroupIndex; 
            }
            else if( m_vecCandidateSpaceGroups[nBestCentricSpaceGroupIndex].fGetProbability() == 
                     m_vecCandidateSpaceGroups[nBestAcentricSpaceGroupIndex].fGetProbability() )
            {
                if(  m_vecCandidateSpaceGroups[nBestCentricSpaceGroupIndex].nGetIndexInLaueClassTable() > 
                     m_vecCandidateSpaceGroups[nBestAcentricSpaceGroupIndex].nGetIndexInLaueClassTable() )
                {
                    m_nSelectedCandidateSpaceGroupIndex = nBestCentricSpaceGroupIndex; 
                }
            }
        }
    }
    else
        return; //error

    //mrpbegin
    if(nFavoredSpaceGroupNumber > 0 && m_vecCandidateSpaceGroups.size() > 1)
    {
        for(size_t i=0; i<m_vecCandidateSpaceGroups.size(); i++)
        {
            if(m_vecCandidateSpaceGroups[i].nGetIntTablesNumber() == nFavoredSpaceGroupNumber)
            {
                m_nSelectedCandidateSpaceGroupIndex = i;
                break;
            }
        }
    }
    //mrpend
}  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find the best candidate space group based on:
// - User requirement for the space group be centric/acentric
// - "Probability" of the space group, according to the Cambridge Structural Database.
// - Index of the corresponding space group entry in the Laue class table
int Cspacegroup_check::nFindBestCandidateSpaceGroup(SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eType)
{
    double      fMaxProbability = 0.0;
    int         nBest = -1;
    
    for(int jj=0; jj < (int)m_vecCandidateSpaceGroups.size(); jj++)
    {
        if( eType == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric &&
            !m_vecCandidateSpaceGroups[jj].bIsCentric() )
        {
            continue;
        }
        else if( eType == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric &&
            m_vecCandidateSpaceGroups[jj].bIsCentric() )
        {
            continue;
        }
        
        if( m_vecCandidateSpaceGroups[jj].fGetProbability() > fMaxProbability )
        {
            fMaxProbability = m_vecCandidateSpaceGroups[jj].fGetProbability();
            nBest = jj;
        }
        else if( nBest > -1 &&
                 m_vecCandidateSpaceGroups[jj].fGetProbability() == fMaxProbability && 
                 m_vecCandidateSpaceGroups[jj].nGetIndexInLaueClassTable() >
                 m_vecCandidateSpaceGroups[nBest].nGetIndexInLaueClassTable() )
        {
            nBest = jj;
        }
    }
    
    return nBest;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vFindProbabilitiesForCandidateSpaceGroups()
{
    // Find the total frequency of selected spacegroups and the spacegroup with the highest frequency
    double      fTotFreq = 0.0;
    int         jj = 0;
    for(jj=0; jj < (int)m_vecCandidateSpaceGroups.size(); jj++)
    {
        fTotFreq +=(double)m_vecCandidateSpaceGroups[jj].nGetCSDBFrequency();
    }
    
    // Set probability to each selected spacegroup 
    for(jj=0; jj < (int)m_vecCandidateSpaceGroups.size(); jj++)
    {
        m_vecCandidateSpaceGroups[jj].vSetProbability((double)m_vecCandidateSpaceGroups[jj].nGetCSDBFrequency() /
                                                      fTotFreq * 100.0);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vMarkInconclusiveConditionsForCandidateSpaceGroupTableEntries
                       (std::vector<CSpaceGroupTableEntry*>& vecCandidates)
{
    for(int ii=0; ii < (int)vecCandidates.size(); ii++)
    {
        if( !vecCandidates[ii]->bIsSelected() )
            continue;

        vecCandidates[ii]->vMarkInconclusiveCondtions(m_vecReflnConditionPtrs);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vSelectRestrictionsFromSelectedCandidateSpaceGroup()
{
    for(int ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
        m_vecReflnConditionPtrs[ii]->vUnselectAllRestrictions();

    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return;
    
    int     nSel = m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].nGetIndexInLaueClassTable();

    m_vecSpaceGroupTableEntryPtrs[nSel]->vSelectRestrictionsInConditions(m_vecReflnConditionPtrs);       
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vGetFullListRestrictionNameByIndex(int iIndex, Cstring& sRestrName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    
    poInfo->bGetRestrictionNameByIndex(iIndex, sRestrName);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vGetFullListRestrictionLongNameByName(Cstring& sRestrName, Cstring& sRestrLongName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);

    poInfo->bGetRestrictionLongNameByName(sRestrName, sRestrLongName);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vGetFullListRestrictionInvertedLongNameByName(Cstring& sRestrName, Cstring& sRestrInvertedLongName)
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    
    poInfo->bGetRestrictionInvertedLongNameByName(sRestrName, sRestrInvertedLongName);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetTotalNumberOfFullListRestrictions()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    
    return poInfo->nGetTotalNumberOfFullListRestrictions();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring Cspacegroup_check::sGetLaueName()
{
    CSpaceGroupTableInfo*   poInfo = DoSpacegroupTableInfo(DTREK_DSGTI_GET);
    
    return poInfo->sGetLaueName(m_stCtrlInfo.eLaue);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return the number of generated solutions
int Cspacegroup_check::nGetSpaceGroupSolutionsFromSpaceGroupTableEntry(int iIndex, 
                                                                       SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eCent)
{
    bool    bCentric    = eCent == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric ? true : false;
    bool    bChiralOnly = m_stCtrlInfo.eChirality == SPACEGROUP_CHECK_CTRL_INFO::enChirality_Chiral ? true : false;
    
    return m_vecSpaceGroupTableEntryPtrs[iIndex]->nAddSpaceGroupSolutionsToVector(m_vecCandidateSpaceGroups,
                                                                                  m_stCtrlInfo.eLaue,
                                                                                  bCentric,
                                                                                  bChiralOnly);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vPrintSystematicAbsencesForSelectedSpaceGroup(Cstring& sOut)
{
    sOut = "";  // just in case
    
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return;

    Cstring     sTemp("");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // First of all we need to establish the maximum length of the "HKL restriction" string for each HKL condition.
    // That way we will know how wide the table will be.
    std::vector<CReflectionRestriction*>    vecFirstSelectedRestrictionPtrs;
    CReflectionRestriction*                 poRestr = NULL;
    
    // Create and initialize a vector to store the maximum length for each condition.
    std::vector<int>                        anMaximumInvertedRestrictionLongNameLength;
    int     ii = 0;
    for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        anMaximumInvertedRestrictionLongNameLength.push_back(0);
    }
   
    // Select all restrictions to be printed out.
    vSelectRestrictionsFromSelectedCandidateSpaceGroup();
    
    while( 0 < nGetTotalNumberOfSelectedRestrictionsInAllConditions() )
    {
        vecFirstSelectedRestrictionPtrs.clear();
        for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
        {
            poRestr = m_vecReflnConditionPtrs[ii]->poGetFirstSelectedRestrictionPtr();
            vecFirstSelectedRestrictionPtrs.push_back(poRestr);
        }

        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            if( NULL != vecFirstSelectedRestrictionPtrs[ii] )
            {
                vecFirstSelectedRestrictionPtrs[ii]->vGetInvertedLongName(sTemp);
                if( (int)sTemp.GetLength() > anMaximumInvertedRestrictionLongNameLength[ii] )
                    anMaximumInvertedRestrictionLongNameLength[ii] = (int)sTemp.GetLength();
                
                vecFirstSelectedRestrictionPtrs[ii]->vSelect(false);
            }
        }
    }
    
    // Restore selections to get ready for the actual printout
    vSelectRestrictionsFromSelectedCandidateSpaceGroup();
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////
    // Print the table header
    int nSGN = nGetSelectedCandidateSpaceGroupIntTablesNumber();
    printf("\nSystematic absences for spacegroup #%d %s (<I/sig(I)> <= %4.2f)\n", 
           nSGN, 
           sGetSelectedCandidateSpaceGroupPresentation().string(),
           m_stCtrlInfo.fIOverSigmaTol);
    /////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Figure out how long the underscoring line should be and print it
    Cstring     sUnderscore("");
    int         nTotalRowLength = 13;  // 13 is the number of positions for the row caption like "Conditions   "
    for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        nTotalRowLength += max(anMaximumInvertedRestrictionLongNameLength[ii], 8);
    }
    nTotalRowLength += (int)m_vecReflnConditionPtrs.size() - 1; // to account for the spaces between the names in a row

    for(ii=0; ii < nTotalRowLength; ii++)
        sUnderscore += '-';
        
    printf("%s\n", sUnderscore.string());
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    char        pcBuf[80];

    Cstring     sTextFormat    ("s ");
    Cstring     sIntegerFormat ("d ");
    Cstring     sFloatFormat   (".2f ");

    Cstring     sFormat    ("");

    printf("Conditions   ");
    
    sOut += " Refln._Type ";
    
    for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        sTemp = m_vecReflnConditionPtrs[ii]->sGetID();
        
        sFormat = "%";
        sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
        sFormat += sTextFormat;
        
        //sprintf(pcBuf, "%8s ", sTemp.string());
        sprintf(pcBuf, sFormat.string(), sTemp.string());

        sOut += pcBuf;
        printf(pcBuf);
    }
    
    
    printf("\n");
    printf("%s\n", sUnderscore.string());
    printf("Total Reflns ");
    sOut += " \n Total_Reflns ";
    
    for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        sFormat = "%";
        sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
        sFormat += sIntegerFormat;
        
        //sprintf(pcBuf, "%8d ", m_vecReflnConditionPtrs[ii]->nGetNumberOfReflns());
        sprintf(pcBuf, sFormat.string(), m_vecReflnConditionPtrs[ii]->nGetNumberOfReflns());
        
        sOut += pcBuf;
        printf(pcBuf);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    while( 0 < nGetTotalNumberOfSelectedRestrictionsInAllConditions() )
    {
        vecFirstSelectedRestrictionPtrs.clear();
        
        for(ii = 0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
        {
            poRestr = m_vecReflnConditionPtrs[ii]->poGetFirstSelectedRestrictionPtr();
            vecFirstSelectedRestrictionPtrs.push_back(poRestr);
        }
    
        //////////////////////////////////////////////////////////////////////////////
        //printf("\n\n  Restriction ");
        printf("\n\nRestrictions ");
        sOut += " \n Restriction ";
        
        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            sFormat = "%";
            sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
            sFormat += sTextFormat;
            
            //if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
                //sprintf(pcBuf,"%8s ","----");
            //else
            //{
                //vecFirstSelectedRestrictionPtrs[ii]->vGetInvertedLongName(sTemp);
                //sprintf(pcBuf,"%8s ", sTemp.string());
            //}

            if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
                sprintf(pcBuf, sFormat.string(),"----");
            else
            {
                vecFirstSelectedRestrictionPtrs[ii]->vGetInvertedLongName(sTemp);
                sprintf(pcBuf, sFormat.string(), sTemp.string());
            }
            
            printf(pcBuf);
            sOut += pcBuf;
        }
        /////////////////////////////////////////////////////////////////////////////
        printf("\nNum Reflns   ");
        sOut += " \n Num_In_Restriction ";
        
        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            
            //if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
                //sprintf(pcBuf,"%8s ","----");
            //else
                //sprintf(pcBuf,"%8d ", vecFirstSelectedRestrictionPtrs[ii]->nGetTotalNumberOfNonCompliantReflns() );
 
            if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sTextFormat;
                sprintf(pcBuf, sFormat.string(), "----");
            }
            else
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sIntegerFormat;
                sprintf(pcBuf, sFormat.string(), vecFirstSelectedRestrictionPtrs[ii]->nGetTotalNumberOfNonCompliantReflns() );
            }

            printf(pcBuf);
            sOut += pcBuf;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\nNum Observed ");
        sOut += " \n Num_I/sig_>=_";
        sprintf(pcBuf, "%.1f", m_stCtrlInfo.fIOverSigmaTol);
        sOut += pcBuf;
        sOut += " ";
        
        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            //if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
                //sprintf(pcBuf,"%8s ","----");
            //else
                //sprintf(pcBuf,"%8d ", vecFirstSelectedRestrictionPtrs[ii]->nGetNumberOfNonCompliantIntegratedSignificantReflns());

            if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sTextFormat;
            
                sprintf(pcBuf, sFormat.string(), "----");
            }
            else
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sIntegerFormat;
            
                sprintf(pcBuf, sFormat.string(), 
                               vecFirstSelectedRestrictionPtrs[ii]->nGetNumberOfNonCompliantIntegratedSignificantReflns());
            }

            printf(pcBuf);
            sOut += pcBuf;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //printf("\n   <I/sig(I)> ");
        printf("\n<I/sig(I)>   ");
        sOut += " \n <I/sig(I)> ";
        
        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            //if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
                //sprintf(pcBuf,"%8s ","----");
            //else
                //sprintf(pcBuf,"%8.2f ", vecFirstSelectedRestrictionPtrs[ii]->fGetIoverSigmaOfNonCompliantIntegratedReflns());

            if( NULL == vecFirstSelectedRestrictionPtrs[ii] )
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sTextFormat;
                sprintf(pcBuf, sFormat.string(), "----");
            }
            else
            {
                sFormat = "%";
                sFormat += anMaximumInvertedRestrictionLongNameLength[ii] <= 8 ? 8 : anMaximumInvertedRestrictionLongNameLength[ii];
                sFormat += sFloatFormat;
                sprintf(pcBuf, sFormat.string(), 
                               vecFirstSelectedRestrictionPtrs[ii]->fGetIoverSigmaOfNonCompliantIntegratedReflns());
            }
            
            printf(pcBuf);
            sOut += pcBuf;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        sOut += " \n ";  
        
        for(ii = 0; ii < (int)vecFirstSelectedRestrictionPtrs.size(); ii++)
        {
            if( vecFirstSelectedRestrictionPtrs[ii] != NULL )
                vecFirstSelectedRestrictionPtrs[ii]->vSelect(false);
        }
    }
    
    printf("\n");
    printf("%s\n", sUnderscore.string());
    printf("\n\n");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vPrintAllConditionsRestrictions(Cstring& sOut)
{
    sOut = "";
    int     ii = 0;
    
    Cstring     sTemp("");
    char        pcBuf[80];

    printf("\nALL Systematic Absences\n");
    
    // We need to establish the length of a row and generate an underscoring string
    int     nRowLength = 14 + (int)m_vecReflnConditionPtrs.size() * 9 - 1; // 14 positions is for row name 
    Cstring     sUnderscore("");
    for(ii=0; ii < nRowLength; ii++)
        sUnderscore += '-';
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("%s\n", sUnderscore.string());
    printf(" Refln. Type  ");
    sOut+= " Refln._Type  ";
    
    int     jj = 0;
    int     kk = 0;
    for(ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        sTemp = m_vecReflnConditionPtrs[ii]->sGetID();
        sprintf(pcBuf,"%8s ",sTemp.string());

        sOut += pcBuf;
        printf(pcBuf);
    }
    
    printf("\n");
    printf("%s\n", sUnderscore.string());
    printf(" Total Reflns ");
    sOut += " \n Total_Reflns ";
    
    for (ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        sprintf(pcBuf, "%8d ", m_vecReflnConditionPtrs[ii]->nGetNumberOfReflns());
        
        sOut+= pcBuf;
        printf(pcBuf);
    }
    
    printf("\n\n");
    //////////////////////////////////////////////////////////////////////////////////////////

    CReflectionCondition*       poCond = NULL;
    CReflectionRestriction*     poRestr = NULL;

    DTREK_DWORD       dwID = 0L;
    
    Cstring                 sRestrName("");
    Cstring                 sRestrLongName("");
    Cstring                 sRestrInvertedLongName("");
    char                    cChar = '\0';

    for(ii=0; ii < nGetTotalNumberOfFullListRestrictions(); ii++)
    {
        vGetFullListRestrictionNameByIndex(ii, sRestrName);
        vGetFullListRestrictionLongNameByName(sRestrName, sRestrLongName);
        vGetFullListRestrictionInvertedLongNameByName(sRestrName, sRestrInvertedLongName);
        
        for(jj=0; jj < 6; jj++)
        {
            if( jj == 3 ) 
            {
                sprintf(pcBuf," %12s ", sRestrLongName.string());
                sOut+=pcBuf;
            } 
            else if( jj == 0 )
            {
                sprintf(pcBuf," %12s ", sRestrInvertedLongName.string());
                sOut+=pcBuf;
            } 
            else
            {
                sprintf(pcBuf," %12s ","");
                sOut+=" notitle ";
            }                    
            printf(pcBuf);
            
            for(kk=0; kk < (int)m_vecReflnConditionPtrs.size(); kk++)
            {
                poCond = m_vecReflnConditionPtrs[kk];
                poRestr = poCond->poGetRestrictionPtrByName(sRestrName);
                
                if( NULL == poRestr )
                    sprintf(pcBuf,"%8s ","------");
                else 
                {
                    switch( jj ) 
                    {
                    case 0: 
                        sprintf(pcBuf,"%8d ", poRestr->nGetTotalNumberOfNonCompliantReflns()); 
                        break;
                    case 1: 
                        sprintf(pcBuf,"%8d ", poRestr->nGetNumberOfNonCompliantIntegratedSignificantReflns());
                        break;
                    case 2: 
                        if( poRestr->bIsInconclusive() )
                            cChar = '?';
                        else if( poRestr->bIsFailed() )
                            cChar = ' ';
                        else
                            cChar = '*'; // conclusive absence
                        sprintf(pcBuf,"%7.2f%c ", poRestr->fGetIoverSigmaOfNonCompliantIntegratedReflns(), cChar); 
                        break;
                    case 3: 
                        sprintf(pcBuf,"%8d ", poRestr->nGetTotalNumberOfCompliantReflns()); 
                        break;
                    case 4: 
                        sprintf(pcBuf,"%8d ", poRestr->nGetNumberOfCompliantIntegratedSignificantReflns()); 
                        break;
                    case 5:
                        sprintf(pcBuf,"%7.2f  ", poRestr->fGetIoverSigmaOfCompliantIntegratedReflns()); 
                        break;
                    }
                }
                
                sOut += pcBuf;
                printf(pcBuf);
            }
            
            printf("\n");
            sOut+= " \n ";
            if( !((jj+1) % 3) )
                printf("\n");
        }
    }
    
    printf("%s\n", sUnderscore.string());
    printf("* This group of reflections is systematically absent\n");
    printf("? Cannot determine whether this group of reflections is systematically absent.\n");
    printf("\n\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vPrintSelectedSpaceGroups(Cstring& sOut)
{
    sOut = "";
    
    if( m_vecCandidateSpaceGroups.size() == 0 )
        return; // nothing to print
    
    Cstring     sLaue = sGetLaueName();
    Cstring     sChir = m_stCtrlInfo.eChirality != SPACEGROUP_CHECK_CTRL_INFO::enChirality_Unknown ? 
                        "chiral" : "not specified";

    // Print the header
    printf("Spacegroups found (Laue class: %s, Chirality: %s)\n", sLaue.string(), sChir.string());
    printf("-------------------------------------------------------------\n");
    //printf(" %6s %8s %12s %10s %10s\n","Number","Name","Presentation","Centricity","Frequency");
    printf(" %6s %8s %12s %10s %11s %8s\n","Number","Name","Presentation","Centricity","Probability", "Absences");
    printf("-------------------------------------------------------------\n");
    
    bool    bIsNonStandard = false;
    char    pcBuf[80];

    for(int ii=0; ii < (int)m_vecCandidateSpaceGroups.size(); ii++)
    {
        if( !m_vecCandidateSpaceGroups[ii].bIsStandard() )
            bIsNonStandard = true;
        
        //sprintf(pcBuf, " %6d %8s %12s %10s %8.2f %% \n",
        sprintf(pcBuf, " %6d %8s %12s %10s %11.1f %8s \n",
                       m_vecCandidateSpaceGroups[ii].nGetIntTablesNumber(),
                       m_vecCandidateSpaceGroups[ii].sGetIntTablesName().string(),
                       m_vecCandidateSpaceGroups[ii].sGetPresentation().string(),
                       m_vecCandidateSpaceGroups[ii].bIsCentric() ? "Centric" : "Acentric",
                       m_vecCandidateSpaceGroups[ii].fGetProbability(),
                       sGetAbsenceTestSummary(ii).string());
        sOut += pcBuf;
        
        printf(pcBuf);
    }

    printf("-------------------------------------------------------------\n");
    printf("\n");
    
    // Print the explanation of the format.
    printf("Probability column shows the relative probability (in percent) for each of the\n"
           "found spacegroups based on the frequencies of just these spacegroups in the\n"
           "small molecule Cambridge Structural Database.\n\n");

    printf("Absences column summarizes systematic absences for reflection conditions shown\n"
           "in the \"Systematic Absences\" table. 'A' stands for a conclusive systematic \n"
           "absence of non-restricted reflections in a condition. '?' means\n"
           "inconclusiveness because of too few reflections in the file in a condition.\n"
           "'-' means a condition is not restricted by the spacegroup.\n\n");

    if( bIsNonStandard ) 
        printf("Spacegroups whose 'Presentation' field differs from the 'Name' field were\n"
               "found in a non-standard presentation. Re-indexing will be required to convert\n"
               "to a standard presentation.\n\n");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetTotalNumberOfSelectedRestrictionsInAllConditions()
{
    int nCount = 0;
    for(int ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        nCount += m_vecReflnConditionPtrs[ii]->nGetNumberOfSelectedRestrictions();
    }
    
    return nCount;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetSelectedCandidateSpaceGroupIntTablesNumber()
{
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return -1;
    
    return m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].nGetIntTablesNumber();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SPACEGROUP_CHECK_CTRL_INFO::eCentricityType Cspacegroup_check::eGetSelectedCandidateSpaceGroupCentricity()
{
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Unknown; 
    
    return (m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].bIsCentric() ? 
            SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric :
            SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetFirstCandidateSpaceGroupIntTablesNumber()
{
    if( 0 == m_vecCandidateSpaceGroups.size() )
        return 0; 

    return m_vecCandidateSpaceGroups[0].nGetIntTablesNumber();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetCandidateSpaceGroupIndexBySGNumber(int nSG)
{
    for(int ii=0; ii < (int)m_vecCandidateSpaceGroups.size(); ii++)
    {
        if( nSG == m_vecCandidateSpaceGroups[ii].nGetIntTablesNumber() )
            return ii;
    }
    
    return -1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SPACEGROUP_CHECK_CTRL_INFO::eCentricityType Cspacegroup_check::eGetCandidateSpaceGroupCentricity(int iIndex)
{
    if( iIndex < 0 || iIndex > (int)m_vecCandidateSpaceGroups.size() - 1)
        return SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Unknown; 
    
    return (m_vecCandidateSpaceGroups[iIndex].bIsCentric() ? 
            SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric :
            SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vSetSelectedCandidateSpaceGroupIndexBySGNumber(int nSG)
{
    int iIndex = nGetCandidateSpaceGroupIndexBySGNumber(nSG);

    vSetSelectedCandidateSpaceGroupIndex(iIndex);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring Cspacegroup_check::sGetSelectedCandidateSpaceGroupPresentation()
{
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return Cstring("");

    return m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].sGetPresentation();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring Cspacegroup_check::sGetSelectedCandidateSpaceGroupName()
{
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return Cstring("");

    return m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].sGetName();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetIndexOfSpaceGroupTableEntryContainingSpaceGroupInStandardSetting(int nSG)
{
    for(int ii=0; ii < (int)this->m_vecSpaceGroupTableEntryPtrs.size(); ii++)
    {
        if( m_vecSpaceGroupTableEntryPtrs[ii]->bIsSpaceGroupInStandardSetting(nSG) )
            return ii;
    }
    
    return -1; // error: space group not found. 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cspacegroup_check::vFillAbsensesArrayOldStyleFromSpaceGroupTableEntry(int iEntry, 
                                                                           int nAbsences[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4])
{
    m_vecSpaceGroupTableEntryPtrs[iEntry]->vFillAbsensesArrayOldStyle(m_vecReflnConditionPtrs, nAbsences);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cspacegroup_check::nGetIndexOfSpaceGroupTableEntryBySelectedCandidateIndex()
{
    if( m_nSelectedCandidateSpaceGroupIndex < 0 || m_nSelectedCandidateSpaceGroupIndex > (int)m_vecCandidateSpaceGroups.size() - 1 )
        return -1;

    return m_vecCandidateSpaceGroups[m_nSelectedCandidateSpaceGroupIndex].nGetIndexInLaueClassTable();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring Cspacegroup_check::sGetAbsenceTestSummary(int iCandidateSpaceGroupIndex)
{
    int iSpaceGroupEntryIndex = m_vecCandidateSpaceGroups[iCandidateSpaceGroupIndex].nGetIndexInLaueClassTable();

    return m_vecSpaceGroupTableEntryPtrs[iSpaceGroupEntryIndex]->sGetAbsenceTestSummary();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is necessary when some of the conditions are inconclusive. We have to compare entries
// between themselves to see if one entry has everything another entry has plus at least one conclusive condition.
// If that is the case then the worse entry will be unselected.
void Cspacegroup_check::vDeselectCandidateSpaceGroupTableEntriesByComparison(std::vector<CSpaceGroupTableEntry*>& vecCandidates)
{
    bool    bEntry_ii_better = false;
    bool    bEntry_jj_better = false;

    bool    bUncomparableEntries = false;

    enumTestStatus  eTestStatus_ii = enUnknown;
    enumTestStatus  eTestStatus_jj = enUnknown;

    Cstring     sRestriction_ii("");
    Cstring     sRestriction_jj("");
    
    for(int ii=0; ii < (int)vecCandidates.size(); ii++)
    {
        if( !vecCandidates[ii]->bIsSelected() )
            continue;
        
        for(int jj=0; jj < (int)vecCandidates.size(); jj++)
        {
            if( jj <= ii )
                continue; // because we have already compared these ii and jj entries
            
            if( !vecCandidates[jj]->bIsSelected() )
                continue;
            
            bEntry_ii_better = false;
            bEntry_jj_better = false;
            
            bUncomparableEntries = false;
            
            for(int kk=0; kk < (int)m_vecReflnConditionPtrs.size(); kk++)
            {
                eTestStatus_ii = vecCandidates[ii]->eGetConditionStatus(m_vecReflnConditionPtrs[kk]); 
                eTestStatus_jj = vecCandidates[jj]->eGetConditionStatus(m_vecReflnConditionPtrs[kk]); 
                
                sRestriction_ii = vecCandidates[ii]->sGetConditionRestriction(m_vecReflnConditionPtrs[kk]);
                sRestriction_jj = vecCandidates[jj]->sGetConditionRestriction(m_vecReflnConditionPtrs[kk]);

                if( eTestStatus_ii == enConclusive && eTestStatus_jj == enUnknown )
                    bEntry_ii_better = true;
                else if ( eTestStatus_jj == enConclusive && eTestStatus_ii == enUnknown )
                    bEntry_jj_better = true;
                else if( eTestStatus_ii != enUnknown && eTestStatus_jj != enUnknown && 
                         !(sRestriction_ii == sRestriction_jj)  )
                    bUncomparableEntries = true;
            }
            
            if( !bUncomparableEntries )
            {
                if( bEntry_ii_better && !bEntry_jj_better )
                {
                    vecCandidates[jj]->vSetSelected(false);
                }
                else if( bEntry_jj_better && !bEntry_ii_better )   
                {
                    vecCandidates[ii]->vSetSelected(false);
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
Cstring Cspacegroup_check::sGetListInconclusiveConditionsForCandidateSpaceGroups()
{
    Cstring     sTemp("");
    
    for(int ii=0; ii < (int)m_vecReflnConditionPtrs.size(); ii++)
    {
        if( m_vecReflnConditionPtrs[ii]->bIsInconclusiveForCandidateSpaceGroup() )
        {
            sTemp += m_vecReflnConditionPtrs[ii]->sGetID();
            sTemp += ' ';
        }
    }
    
    return sTemp;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Try to convert one set of absences into another set of absences using a transformation.  
// The transformation assumes that the restrictions are ordered as found in the monoclinic table 
// and the orthorhombic table. This is refactored Thad's old code. I haven't really tested all of it.
eIndexTransform Cspacegroup_check::eFindTranslationPermForSelectedSpaceGroup()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // For monoclinics, we must use one of the cell type conversion matrices.  Parse name to see if it is a monoclinic cell.
    Cstring     sSelSGName = sGetSelectedCandidateSpaceGroupName();
    int         nLength = (int)sSelSGName.GetLength();

    if( sSelSGName.GetAt(nLength-1) == ')' && sSelSGName.GetAt(nLength-3) == 'C' ) 
    {
        return (eIndexTransform)((sSelSGName.GetAt(nLength-2) - '2') + TRANS_C21);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int     iSpaceGroupSelTableEntry = nGetIndexOfSpaceGroupTableEntryBySelectedCandidateIndex();
    int     nSelSGNumberIntTables = nGetSelectedCandidateSpaceGroupIntTablesNumber();
    int     iStandardSpaceGroupTableEntry = nGetIndexOfSpaceGroupTableEntryContainingSpaceGroupInStandardSetting(nSelSGNumberIntTables);
        
    int     nRestrictionsOrig[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4] = {0};
    int     nRestrictionsTarg[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4] = {0};

    vFillAbsensesArrayOldStyleFromSpaceGroupTableEntry(iSpaceGroupSelTableEntry, nRestrictionsOrig); 
    vFillAbsensesArrayOldStyleFromSpaceGroupTableEntry(iStandardSpaceGroupTableEntry, nRestrictionsTarg);

    enCond eCondMap[6][7] =
    {
        {enCond_HKL,   enCond_0KL,    enCond_H0L,    enCond_HK0,    enCond_H00,    enCond_0K0,    enCond_00L},   /*ABC*/ 
        {enCond_HKL,   enCond_H0L,    enCond_HK0,    enCond_0KL,    enCond_0K0,    enCond_00L,    enCond_H00},   /*CAB*/ 
        {enCond_HKL,   enCond_HK0,    enCond_0KL,    enCond_H0L,    enCond_00L,    enCond_H00,    enCond_0K0},   /*BCA*/ 
        {enCond_HKL,   enCond_H0L,    enCond_0KL,    enCond_HK0,    enCond_0K0,    enCond_H00,    enCond_00L},   /*BAC*/ 
        {enCond_HKL,   enCond_0KL,    enCond_HK0,    enCond_H0L,    enCond_H00,    enCond_00L,    enCond_0K0},   /*ACB*/ 
        {enCond_HKL,   enCond_HK0,    enCond_H0L,    enCond_0KL,    enCond_00L,    enCond_0K0,    enCond_H00}    /*CBA*/ 
    };

    enRestr eRestrMap[6][13] =
    {
        {enRestr_H, enRestr_K, enRestr_L, enRestr_KL, enRestr_HL, enRestr_HK, enRestr_4H,  enRestr_4K,  enRestr_4L,  enRestr_4KL, enRestr_4HL, enRestr_4HK,    enRestr_HKL},        /*ABC*/ 
        {enRestr_K, enRestr_L, enRestr_H, enRestr_HL, enRestr_HK, enRestr_KL, enRestr_4K,  enRestr_4L,  enRestr_4H,  enRestr_4HL, enRestr_4HK, enRestr_4KL,    enRestr_HKL},        /*CAB*/ 
        {enRestr_L, enRestr_H, enRestr_K, enRestr_HK, enRestr_KL, enRestr_HL, enRestr_4L,  enRestr_4H,  enRestr_4K,  enRestr_4HK, enRestr_4KL, enRestr_4HL,    enRestr_HKL},        /*BCA*/ 
        {enRestr_K, enRestr_H, enRestr_L, enRestr_HL, enRestr_KL, enRestr_HK, enRestr_4K,  enRestr_4H,  enRestr_4L,  enRestr_4HL, enRestr_4KL, enRestr_4HK,    enRestr_HKL},        /*BAC*/ 
        {enRestr_H, enRestr_L, enRestr_K, enRestr_KL, enRestr_HK, enRestr_HL, enRestr_4H,  enRestr_4L,  enRestr_4K,  enRestr_4KL, enRestr_4HK, enRestr_4HL,    enRestr_HKL},        /*ACB*/ 
        {enRestr_L, enRestr_K, enRestr_H, enRestr_HK, enRestr_HL, enRestr_KL, enRestr_4L,  enRestr_4K,  enRestr_4H,  enRestr_4HK, enRestr_4HL, enRestr_4KL,    enRestr_HKL}         /*CBA*/ 
    };

    int     nx = 0;
    int     ny = 0;

    int     nTrans = 0;     // Transformation to be used
    int     nCond0 = 0;     // Index in nRestrictions array for conditions as they are compared
    int     nCond1 = 0;      
    int     nCondCt0 = 0;   // Count of Restrictions
    int     nCondCt1 = 0;    

    bool bWorks = false;
    for(nTrans = 0; nTrans < 6; nTrans++)
    {
        for(nCond0 = 0; nCond0 < 7; nCond0++)
        {
            // The condition after translation.
            nCond1=(int) eCondMap[nTrans][nCond0];
            
            // Discover the number of Restrictions in the condition for the source and destination.
            for(nCondCt0 = 0; nCondCt0 < 4; nCondCt0++) 
                if( nRestrictionsOrig[nCond0][nCondCt0] == enRestr_Unknown )
                    break;
            
            for(nCondCt1=0; nCondCt1 < 4; nCondCt1++) 
            {
                if( nRestrictionsTarg[nCond1][nCondCt1] == enRestr_Unknown )
                    break;
            }

            // The number of Restrictions in the source and destination differ.
            if( nCondCt1 != nCondCt0 )
                break;
            
            // The number of Restrictions are the same.  See if the Restrictions are equivalent.
            // (They might be in a differnt order)
            for(nx = 0, bWorks = true; nx < nCondCt0; nx++)
            {
                for(ny = 0; ny < nCondCt1; ny++)
                {
                    enRestr  eVar =(enRestr)nRestrictionsTarg[nCond1][ny];
                    
                    if( eRestrMap[nTrans][nRestrictionsOrig[nCond0][nx]] == eVar )
                        break;
                }

                // Did we get all the way through the list of Restrictions for THIS Restriction?
                if( ny == nCondCt1 )
                {
                    bWorks = false;
                    break;
                }
            }
            
            // Did we get all the way through the list of Restrictions for all Restrictions?
            if( !bWorks ) 
                break;
        }
        
        // Did we get all the way through the list of conditions?
        if( nCond0 == 7 ) 
            break;
    }

    if( nTrans == 6 ) 
    {
        REPORT_ERROR_FILE_LINE;
    }

    return (eIndexTransform)nTrans;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
