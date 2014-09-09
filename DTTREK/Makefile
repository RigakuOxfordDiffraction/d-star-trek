# This is a Makefile for gcc
# 01-Mar-95	JWP	Molecular Structure Corporation
#
#
# Copyright (C) 2014 Rigaku Americas Corporation
#                    9009 New Trails Drive
#                    The Woodlands, TX, USA  77381
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice(s), this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice(s), this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the Rigaku Americas Corporation nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#

# 
# Define method to get system architecture to use in subsequent file names:
#

SHELL	   =	/bin/sh
GPP_DIR    =	/usr/freeware/bin/
#GPP_DIR    =	/usr/local/bin/
CC	   =	$(GPP_DIR)g++ -DIRIX -Wno-non-template-friend -DBINARY_NOT_BIN
#On linux anyways, -O2 breaks dtsvd.o
OPTLEV     =    -g -O1
#OPTLEV     =    -g
#GINCLUDEDIR=	/usr/local/lib/g++-include
#GINCLUDEDIR=	/usr/local/include/g++-3
GINCLUDEDIR=	.

#############################################################################

# For RAXIS     support add -DRAXIS     to CCDEFS and raxis.o  to MSCOBJS below.
# For SIEMENS   support add -DSIEMENS   to CCDEFS and bs.o     to MSCOBJS below.
# For MARCCD    support add -DMARCCD    to CCDEFS and marccd.o to MSCOBJS below.
# For MARIP     support add -DMARIP     to CCDEFS and marccd.o to MSCOBJS below.
# For ADSCCCD   support add -DADSCCCD   to CCDEFS and adscccd.o to MSCOBJS below
# For BrandeisCCD .port add -DBRANDEISCCD  CCDEFS and brandeisccd.o MSCOBJS bel
# For DIP2030   support add -DDIP2030   to CCDEFS and marccd.o to MSCOBJS below.
# For PARALLEL  support add -DPARALLEL  to CCDEFS and set ANL_DTREKIO_INC below.
# For ANL_TIMER support add -DANL_TIMER to CCDEFS and set ANL_TIMER_INC   below.

DZ           =  -DDEMO_LICENSE
#DZ           =
CCDEFS       = $(DZ)

# Choose only ONE of the next definitions

MSCOBJSDZ   = des2.o raxis.o bs.o marccd.o dip2030.o adscccd.o medoptics.o brandeisccd.o marip.o pck.o winbmp.o

MSCOBJS    = $(MSCOBJSDZ)

# If -DPARALLEL, you need par_dtrek.h, also set ANL_LDFLAGS and ANL_LIBS below.

ANL_DTREKIO_INC =
#ANL_DTREKIO_INC = par_dtrek.h

# If -DANL_TIMER, you need anlTimer.h, also set ANL_LDFLAGS and ANL_LIBS below.

ANL_TIMER_INC =
#ANL_TIMER_INC = anlTimer.h

# Choose only ONE pair of the next pairs of definitions
# You need the first pair if you have -DANL_TIMER and -DPARALLEL

#ANL_LDFLAGS = -lanl_timer$(OSTYPE)
#ANL_LIBS    = libanl_timer$(OSTYPE).a
#ANL_LDFLAGS = -lanl_timer$(OSTYPE) -lpar_dtrek$(OSTYPE) 
#ANL_LIBS    = libanl_timer$(OSTYPE).a libpar_dtrek$(OSTYPE).a 
#ANL_LDFLAGS = -lpar_dtrek$(OSTYPE) 
#ANL_LIBS    =  libpar_dtrek$(OSTYPE).a 
ANL_LDFLAGS =
ANL_LIBS    =   

#############################################################################

CCFLAGS	   =	$(OPTLEV) $(CCDEFS) -I. -I../DTDISPLAY -I$(GINCLUDEDIR) \
		-L../DTDISPLAY -L.  -L/usr/lib32/mips3

CCFLAGSDBG =	-g $(CCDEFS) -I. -I../DTDISPLAY -I$(GINCLUDEDIR) \
		-L../DTDISPLAY -L.

GLIBDIR	   =	-L/usr/local/lib
COMPILE.cc =	$(CC) -c $(CCFLAGS)
INSTALLDIR =    ../../bin$(DTREK_BINSUFFIX)


ALL	   =	dtfind dtindex dtrefine dtpredict dtstrategy \
                dtintegrate dtprofit dtcorrect dtreflnmerge \
		dtscaleaverage dtorient dtheaderedit \
		dtnonunfedit dtaverage dtrek2raxis dtextractheader \
		dtdezinger dtdarksubtract dtraxis dtcreate dtadd2 dtcompare \
		dtcompose dtdecompose dtbin dtheadermerge \
		dtcell dtsmooth dtdistor \
		dtx2refln dtsbcdark dtsbcdarkmask dtrebatch dtdisplayupdate \
		licenseinfo 

# OSF/Motif is NOT needed to compile modules in this directory
# It helps with distribution if -static is used for linux with these modules

#LIBS	   =	$(MSCOBJS) -ldtrekX$(OSTYPE) -lXt -lX11 -lm
###linux SLIBS      =	-static -lXp -lXext  -lSM -lICE -lglib
LIBS	   = $(ANL_LDFLAGS) -ldtrekX$(OSTYPE) -lXt -lX11 # -lm -lcrypt

LIBOBJS	   = 	Cimage.o Cimage_header.o dskio.o Cstring.o Cscan.o \
		Cgoniometer.o Cnonunf.o BitmapRLE.o Crotation.o Cwavelength.o \
		dtrekvec.o Creflnlist.o Crefln.o Cspacegroup.o Cspatial.o \
		Ccrystal.o  dtreksys.o Cstat.o dtsvd.o dtrekvec.o \
		Csource.o Cdetector.o C3Ddata.o Cgenetic.o CEmpericalProfit.o \
                dtreklogic.o dtarray.o \
                DevErrCode.o CinterpMesh.o CScanBitmap.o Cprofit2D.o $(MSCOBJS)

ANL_PAROBJS =	par_dtrek.o

ANL_TIMOBJS =	anlTimer.o

RM	   = rm -f

all:		$(ALL)

tCimage:	tCimage.o Cimage.o Cimage_header.o Cstring.o Crefln.o \
		Creflnlist.o Ccrystal.o Cspacegroup.o dskio.o \
		dtreklogic.o dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

lib:            libdtrek$(OSTYPE).a  $(ANL_LIBS)

libdtrek$(OSTYPE).a:	$(LIBOBJS)
		$(AR) -rv $@ $^

libpar_dtrek$(OSTYPE).a:	$(ANL_PAROBJS)
		$(AR) -rv $@ $^

libanl_timer$(OSTYPE).a:	$(ANL_TIMOBJS)
		$(AR) -rv $@ $^

main.o:		main.cc Cimage.h Cimage_header.h

Cimage_header.o: Cimage_header.cc Cimage_header.h

Cimage.o:	Cimage.cc  Cimage.h

par_dtrek.o:	par_dtrek.cc par_dtrek.h

Cmem.o:		Cmem.cc Cmem.h

CEmpericalProfit.o:	CEmpericalProfit.cc CEmpericalProfit.h

Cgenetic.o:	Cgenetic.cc Cgenetic.h

Cspatial.o:	Cspatial.cc  Cspatial.h

Cspacegroup.o:	Cspacegroup.cc  Cspacegroup.h Dtrek.h  Cspacegroup_table.h

Ccrystal.o:	Ccrystal.cc  Ccrystal.h Dtrek.h Cimage_header.h

Cgoniometer.o:	Cgoniometer.cc  Cgoniometer.h Dtrek.h Cimage_header.h

Csource.o:	Csource.cc  Csource.h Dtrek.h Cimage_header.h Cwavelength.h

Cdetector.o:	Cdetector.cc  Cdetector.h Dtrek.h Cimage_header.h

Cwavelength.o:	Cwavelength.cc  Cwavelength.h Dtrek.h Cimage_header.h

tCspatial:	tCspatial.o Cspatial.o Cimage.o Cimage_header.o Cstring.o \
		dskio.o Crefln.o Creflnlist.o Cspacegroup.o Csource.o \
		Ccrystal.o Cdetector.o Cnonunf.o BitmapRLE.o Cgoniometer.o Cwavelength.o \
		C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dtsvd.o dtarray.o dtrekvec.o dtreklogic.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCreflnlist:	tCreflnlist.o Creflnlist.o Crefln.o Cimage.o Cimage_header.o \
		Ccrystal.o Cspacegroup.o Cstring.o dtrekvec.o dskio.o \
		dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tR:		tR.o Creflnlist.o Crefln.o Cimage.o Cimage_header.o dskio.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tR.o:		tR.cc Creflnlist.h Crefln.h Dtrek.h

tCscan:		tCscan.o Cscan.o Crotation.o Cnonunf.o BitmapRLE.o Cimage.o \
		Cimage_header.o Cstring.o dtrekvec.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o \
		Ccrystal.o dtreklogic.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o dtreksys.o Cdetector.o \
		Cgoniometer.o Cspatial.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCscan.o:	tCscan.cc Cscan.h

tCdetector:	tCdetector.o Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Cimage.o \
		Cgoniometer.o Csource.o Cwavelength.o Cimage_header.o \
		Cimage.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crefln.o Creflnlist.o Ccrystal.o \
		Cspacegroup.o dtreklogic.o  dtreksys.o \
		Cstring.o dskio.o dtrekvec.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCdetector.o:	tCdetector.cc Dtrek.h Cdetector.h Csource.h Cstring.h

tCnonunf:	tCnonunf.o Cnonunf.o BitmapRLE.o Cimage.o Cimage_header.o dskio.o \
		Cstring.o Cimage.o Cimage_header.o Crotation.o \
		Cspacegroup.o dtreklogic.o dtrekvec.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crefln.o \
		Creflnlist.o Ccrystal.o Cdetector.o Cspatial.o Cgoniometer.o \
		$(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCnonunf.o:	tCnonunf.cc Cnonunf.h 

Cnonunf.o:	Cnonunf.cc Cnonunf.h BitmapRLE.h

BitmapRLE.o:	BitmapRLE.cc BitmapRLE.h

Cnonunf.h:	Cimage.h C3Ddata.h Cstring.h

tCstring.o:	tCstring.cc Cstring.h

tCstring:	tCstring.o Cstring.o 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

Cjim.o:		Cjim.cc Cjim.h

testjim:	testjim.o Cjim.o Cstring.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

Cnonunf.h:	Dtrek.h

Cimage.h:	Cimage_header.h

Cscalemerge.o:	Cscalemerge.cc Cscalemerge.h Creflnlist.h Crefln.h Dtrek.h \
		Cimage_header.h Ccrystal.h

Creflnlist.o:	Creflnlist.cc Creflnlist.h Crefln.h Dtrek.h
		$(CC) -c $(CCFLAGS) Creflnlist.cc

Crefln.o:	Crefln.cc Creflnlist.h Crefln.h Dtrek.h
		$(CC) -c $(CCFLAGS) Crefln.cc

Cinterp.o:	Cinterp.cc

CinterpMesh.o:	CinterpMesh.cc CinterpMesh.h

Claue.o:	Claue.cc dtcell.h

Cabsence.o:	Cabsence.cc dtcell.h

dtcell.o:	dtcell.cc dtcell.h
		$(CC) -c $(CCFLAGS) -DCCELL_DEFINE_DATA dtcell.cc

Cscan.o:	Cscan.cc Cscan.h

Cscan.h:	Crotation.h Cnonunf.h Cspatial.h Cimage.h

Crotation.o:	Crotation.cc Crotation.h Dtrek.h

Crotation.h:	Cstring.h Cimage_header.h dtrekvec.h

anlTimer.o:     anlTimer.cc $(ANL_TIMER_INC)

Cprofit2D.o:	Cprofit2D.cc Cprofit2D.h
		$(CC) -c -o $@ $(CCFLAGS) $<

CScanBitmap.o:	CScanBitmap.cc CScanBitmap.h

Cintegrate.o:	Cintegrate.cc Cintegrate.h Cfind.h Crefine.h Cpredict.h
		$(CC) -DLESSMEM -c -o $@ $(CCFLAGS) $<

Cintegrate.h:	Cscan.h Crefine.h Cpredict.h C3Ddata.h Crefln.h

Cspatial.h:	Cimage.h

Creflnlist.h:	Crefln.h

tCimage.o:	tCimage.cc Cimage.h

dtcell:		dtcell.o Claue.o Cindex.o Cabsence.o Cimage_header.o CCellReduce.o \
		Cstring.o Crefln.o Creflnlist.o Cspacegroup.o dskio.o \
		dtreklogic.o Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o \
		Cwavelength.o Cgoniometer.o Ccrystal.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
		Cspatial.o Cimage.o Crotation.o Cstat.o dtsvd.o dtarray.o \
		Ccrystal.o dtrekvec.o  Cscan.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtheaderedit:	dtheaderedit.o Cimage.o Cimage_header.o Cstring.o Crotation.o \
		 Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		Cgoniometer.o Ccrystal.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtheaderedit.o:	dtheaderedit.cc Cimage.h

dtextractheader: dtextractheader.o Cimage.o Cimage_header.o Cstring.o \
		 Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o Crotation.o Cscan.o \
		Crotation.o Cscan.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtextractheader.o:	dtextractheader.cc Cimage.h

dtcompose:	dtcompose.o Cimage.o Cimage_header.o Cstring.o Crotation.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o Cscan.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtcompose.o:	dtcompose.cc Cimage.h

dtdecompose:	dtdecompose.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS) 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtdecompose.o:	dtdecompose.cc Cimage.h

dtcreate:	dtcreate.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o  dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtcreate.o:	dtcreate.cc Cimage.h

dtbin:		dtbin.o Cimage.o Cimage_header.o Cstring.o Crefln.o \
		Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Csource.o Cwavelength.o Cstat.o \
		Crotation.o Cscan.o \
		Cdetector.o Cgoniometer.o Cspatial.o Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS) 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtbin.o:	dtbin.cc Cimage.h

dtdezinger:	dtdezinger.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtdezinger.o:	dtdezinger.cc Cimage.h

dtsbcdarkmask:	dtsbcdarkmask.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^ $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtsbcdark:	dtsbcdark.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o Cstat.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cscan.o Crotation.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^ $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtdarksubtract:	dtdarksubtract.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^ $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtdarksubtract.o: dtbin.cc Cimage.h

dtmakemarks:	dtmakemarks.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Ccrystal.o dtrekvec.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		Crotation.o Cscan.o \
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtmakemarks.o:	dtmakemarks.cc Cimage.h

dtreference:	dtreference.o Cimage.o Cimage_header.o Cstring.o Cspatial.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o Cscan.o Cstat.o \
		Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crotation.o \
		Crotation.o Cscan.o \
                Cnonunf.o BitmapRLE.o  dtreksys.o dtrekvec.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtreference.o:	Cimage.h

dtx2refln.o:	dtx2refln.cc Cstring.h des.h

dtx2refln:	dtx2refln.o Cstring.o dtreksys.o Cspacegroup.o Crefln.o \
 		Creflnlist.o dtreklogic.o dtrekvec.o Cimage_header.o Cimage.o \
		dskio.o Ccrystal.o Cgoniometer.o Cindex.o Creflnlist.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o Crefln.o Cscan.o \
		Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crotation.o CCellReduce.o \
                Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Cstat.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtscalepack2refln.o:	dtscalepack2refln.cc Cstring.h des.h

dtscalepack2refln:	dtscalepack2refln.o Cstring.o dtreksys.o \
		Cspacegroup.o Crefln.o \
 		Creflnlist.o dtreklogic.o dtrekvec.o Cimage_header.o Cimage.o \
		dskio.o Ccrystal.o Cgoniometer.o Cindex.o Creflnlist.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o Crefln.o Cscan.o \
		Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crotation.o \
                Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Cstat.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrefln2raxis:	dtrefln2raxis.o Crefln.o Creflnlist.o Cimage.o
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o dskio.o \
		dtreklogic.o  dtreksys.o dtrekvec.o $(MSCOBJS) 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrefln2raxis.o:	dtrefln2raxis.cc Crefln.h Creflnlist.h Ccrystal.h \
			Cimage_header.h

dtrebatch:	dtrebatch.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o dskio.o \
                Cwavelength.o Cgoniometer.o Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o \
                Csource.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		dtreklogic.o  dtreksys.o dtrekvec.o dtsvd.o dtarray.o $(MSCOBJS) 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrebatch.o:	dtrebatch.cc Crefln.h Creflnlist.h Ccrystal.h \
			Cimage_header.h

dtrefln2numabs:	dtrefln2numabs.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o dskio.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		dtreklogic.o  dtreksys.o dtrekvec.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrefln2numabs.o:	dtrefln2numabs.cc Crefln.h Creflnlist.h Ccrystal.h \
			Cimage_header.h

dtraxis2refln:	dtraxis2refln.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Cspacegroup.o Ccrystal.o Cstring.o dskio.o \
		dtreklogic.o dtrekvec.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtraxis2refln.o:	dtraxis2refln.cc Crefln.h Creflnlist.h Ccrystal.h \
			Cimage_header.h

dtheadermerge:	dtheadermerge.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		Cgoniometer.o Ccrystal.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
		Crotation.o Cscan.o \
		dtreklogic.o dtrekvec.o dskio.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^ $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtheadermerge.o:dtheadermerge.cc Cimage_header.h

dtreflnoverlap:	dtreflnoverlap.o Creflnlist.o Crefln.o \
		Cimage.o Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o \
		Cdetector.o Csource.o Cwavelength.o Cgoniometer.o Cspatial.o \
		Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o dtrekvec.o dskio.o \
		dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtreflnmerge:	dtreflnmerge.o Crefln.o Creflnlist.o Cimage.o Cimage_header.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		Cgoniometer.o Ccrystal.o Cspacegroup.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
                Cstat.o Crotation.o Cscan.o \
		dtsvd.o dtarray.o dtreklogic.o dtrekvec.o dskio.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtreflnmerge.o:	dtreflnmerge.cc Crefln.h Creflnlist.h Ccrystal.h \
		Cspacegroup.h Cimage_header.h

dtscaleaverage.o:	dtscaleaverage.cc Cscaleaverage.h 
		$(CC) -DDISTOR_ONLY -c -o $@ $(CCFLAGS) $<

Cscaleaverage.o:	Cscaleaverage.cc Cscaleaverage.h 
		$(CC) -DDISTOR_ONLY -c -o $@ $(CCFLAGS) $<

Cscale.o:	Cscale.cc Cscaleaverage.h
		$(CC) -DDISTOR_ONLY -c -o $@ $(CCFLAGS) $<

CSphericalHarmonic.o:	CSphericalHarmonic.cc CSphericalHarmonic.h
		$(CC) -c -o $@ $(CCFLAGS) $<

CabsorbModel.o:	CabsorbModel.cc Cabsorb.h Cstring.h Dtrek.h Cstring.h Creflnlist.h
		$(CC) -c -o $@ $(CCFLAGS) $<

CabsorbSurface.o:	CabsorbSurface.cc Cabsorb.h Cstring.h Dtrek.h Cstring.h Creflnlist.h
		$(CC) -c -o $@ $(CCFLAGS) $<

dtscaleaverage:	dtscaleaverage.o Cscaleaverage.o Crefln.o Cscan.o Crotation.o \
		Creflnlist.o Cimage.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		CSphericalHarmonic.o \
		Cgoniometer.o Ccrystal.o Cspacegroup.o \
		dtreklogic.o dtrekvec.o dskio.o dtreksys.o \
		Cscale.o dtsvd.o dtarray.o Cstat.o $(MSCOBJS) 
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtscalemerge:	dtscalemerge.o Cscalemerge.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o Cstat.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		Cgoniometer.o Ccrystal.o Cspacegroup.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
		dtreklogic.o dtrekvec.o dskio.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtscalemerge.o:	dtscalemerge.cc Cscalemerge.h Crefln.h Creflnlist.h \
			Ccrystal.h Cspacegroup.h

dtlin:		dtlin.o Crefln.o Creflnlist.o Cimage.o \
		Cimage_header.o Ccrystal.o Cspacegroup.o Cstring.o \
		dtreklogic.o dtrekvec.o dskio.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtlin.o:	dtlin.cc Crefln.h Creflnlist.h \
			Ccrystal.h Cspacegroup.h

dtindex:	dtindex.o Cindex.o Cscan.o Creflnlist.o Crotation.o \
			CEyeballIndex.o \
			CCellReduce.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtindex.o:	dtindex.cc Dtrek.h Cindex.h dtrekvec.h C3Ddata.h

dtrefine:	dtrefine.o Crefine.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o CScanBitmap.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			Cpredict.o Cfind.o Cstat.o dtsvd.o dtarray.o dtarray.o \
			dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrefine.o:	dtrefine.cc Dtrek.h Crefine.h dtrekvec.h C3Ddata.h

dttransforms:	dttransforms.o Ccalibrate.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cstat.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			Csocket.o CsocketImg.o dtsvd.o dtarray.o \
			dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dttransformc:	dttransformc.o Ccalibrate.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cstat.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			Csocket.o CsocketImg.o dtsvd.o dtarray.o \
			dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

Csocket.o:	Csocket.cc Csocket.h

CsocketImg.o:	CsocketImg.cc CsocketImg.h Csocket.h Cimage.h

dtcalibrate:	dtcalibrate.o Ccalibrate.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cinterp.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtcalibrate.o:	dtcalibrate.cc Dtrek.h Ccalibrate.h dtrekvec.h C3Ddata.h

dtcorrect.o:  dtcalibrate.cc Dtrek.h Ccalibrate.h dtrekvec.h C3Ddata.h
		$(CC) -DCORRECT_ONLY -c -o $@ $(CCFLAGS) $<

Ccorrect.o:  Ccalibrate.cc Dtrek.h Ccalibrate.h dtrekvec.h
		$(CC) -DCORRECT_ONLY -c -o $@ $(CCFLAGS) $<

dtcorrect:	dtcorrect.o Ccorrect.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cinterp.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtdistor:	dtdistor.o Cdistor.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cinterp.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

Cdistor.o:  Ccalibrate.cc Dtrek.h Ccalibrate.h dtrekvec.h
		$(CC) -DDISTOR_ONLY -c -o $@ $(CCFLAGS) $<

dtdistor.o:  dtcalibrate.cc Dtrek.h Ccalibrate.h dtrekvec.h C3Ddata.h
		$(CC) -DDISTOR_ONLY -c -o $@ $(CCFLAGS) $<

dtstrategy:	dtstrategy.o Cstrategy.o Cscan.o Creflnlist.o Crotation.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o \
		Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
		Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
		Cimage_header.o Cstring.o dskio.o dtrekvec.o \
		Cpredict.o Cfind.o Cstat.o dtsvd.o dtarray.o dtarray.o \
		dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtstrategy.o:	dtstrategy.cc Dtrek.h dtrekvec.h Cstrategy.h

dtreduce:	dtreduce.o Cindex.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o dtsvd.o dtarray.o \
			dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtreduce.o:	dtreduce.cc Dtrek.h Cindex.h dtrekvec.h C3Ddata.h

dtnonunf:	dtnonunf.o Cimage.o Cimage_header.o Cstring.o \
		Crefln.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		dtrekvec.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtnonunf.o:	dtnonunf.cc Cimage.h

flatnonunf:	flatnonunf.o Cimage.o Cimage_header.o Cstring.o dskio.o \
		$(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)


tCfind:		tCfind.o Cfind.o Cscan.o Creflnlist.o Crotation.o \
                             Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
                             Cimage_header.o dskio.o Cstring.o \
			     dtrekvec.o Crotation.o Cspacegroup.o Ccrystal.o \
			     dtreklogic.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCfind.o:	tCfind.cc Dtrek.h Cfind.h Cscan.h Creflnlist.h

tCorient:	tCorient.o Corient.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o \
                        Cimage.o  Crefln.o Creflnlist.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCorient.o:	tCorient.cc Dtrek.h Corient.h dtrekvec.h

dtorient:	dtorient.o Corient.o NTHUfixup.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o Crotation.o Cscan.o \
		Ccrystal.o Csource.o Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o Cimage.o  Crefln.o Creflnlist.o \
		Cimage_header.o Cstring.o dskio.o dtrekvec.o \
		dtreklogic.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtorient.o:	dtorient.cc Dtrek.h Corient.h dtrekvec.h

NTHUfixup.o:	NTHUfixup.cc Dtrek.h Cimage_header.h Cstring.h Ccrystal.h

tCpredict:	tCpredict.o Cpredict.o Cscan.o Creflnlist.o Crotation.o \
			    Cspacegroup.o Cwavelength.o Cgoniometer.o \
			    Ccrystal.o Csource.o Cdetector.o \
                            Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
			    Corient.o Cimage_header.o dskio.o dtrekvec.o \
			    dtreklogic.o Cstring.o Cimage.o \
			    Cimage_header.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCpredict.o:	tCpredict.cc Dtrek.h Cpredict.h dtrekvec.h

tCrefine:	tCrefine.o Crefine.o Cscan.o Creflnlist.o Crotation.o \
			Cspacegroup.o Cwavelength.o Cgoniometer.o \
			Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o \
                        Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o \
                        Cimage_header.o Cstring.o dskio.o dtrekvec.o \
			dtreklogic.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCrefine.o:	tCrefine.cc Dtrek.h Crefine.h dtrekvec.h

jtest.o:	jtest.cc Cstring.h des.h
		$(CC) -c -o $@ jtest.cc  -I/usr/include/netman $(CCFLAGSDBG)

jtest:		jtest.o Cstring.o dtreksys.o Cspacegroup.o Crefln.o \
 		Creflnlist.o dtreklogic.o dtrekvec.o Cimage_header.o Cimage.o \
		dskio.o Ccrystal.o Cgoniometer.o Cindex.o Creflnlist.o \
		Cspacegroup.o Cwavelength.o Cgoniometer.o Crefln.o Cscan.o \
		Ccrystal.o Csource.o Cdetector.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Crotation.o \
                Cnonunf.o BitmapRLE.o Cspatial.o Cimage.o Crefln.o Cstat.o dtsvd.o dtarray.o \
		$(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

licenseinfo.o:	licenseinfo.cc
		$(CC) -c -O3 -o $@ licenseinfo.cc  $(CCFLAGS)

licenseinfo:	licenseinfo.o Cstring.o dtreksys.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS) -lm

des2.o:		des2.cc des2.h
		$(CC) -c -O3 -o $@ des2.cc  $(CCFLAGS)

license.o:	license.cc lic1.jj Cstring.h des2.h
		$(CC) -c -O3 -o $@ license.cc  $(CCFLAGS)

license:	license.o des2.o Cstring.o dtreksys.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) -lm -lcrypt

licensekeydemo.o:	license.cc lic1.jj Cstring.h des2.h
		$(CC) -DDEMO_LICENSE -DLICENSE_KEY -c -g -o $@ license.cc  $(CCFLAGS)

licensekeynotdemo.o:	license.cc lic1.jj Cstring.h des2.h
		$(CC)  -DLICENSE_KEY -c -g -o $@ license.cc  $(CCFLAGS)

licensekeynotdemo:	licensekeynotdemo.o des2.o Cstring.o dtreksys.o
		$(CC) -DLICENSE_KEY -o $@ $^  $(CCFLAGS) $(GLIBDIR) -lm -lcrypt

licensekeydemo:	licensekeydemo.o des2.o Cstring.o dtreksys.o
		$(CC) -DLICENSE_KEY -DDEMO_LICENSE -o $@ $^  $(CCFLAGS) $(GLIBDIR) -lm -lcrypt

jjtest:		jtest.o dtrekvec.o Cstring.o dskio.o Cimage.o \
	        Cimage_header.o Crefln.o Creflnlist.o Cspacegroup.o \
		Ccrystal.o  C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o dtreklogic.o dtreksys.o \
		Cscan.o Crotation.o Cgoniometer.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cdetector.o Csource.o Cwavelength.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtplot.o:	dtplot.cc dtrekvec.h Dtrek.h Cstring.h

dtplot:		dtplot.o \
		Creflnlist.o Crotation.o Cspacegroup.o \
                Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cimage.o Crefln.o \
		Cimage_header.o dskio.o Cstring.o dtreklogic.o dtrekvec.o \
		Ccrystal.o Cdetector.o Csource.o Cwavelength.o \
		Cgoniometer.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtfind:		dtfind.o Cfind.o Creflnlist.o Crotation.o Cspacegroup.o \
                Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cimage.o Crefln.o \
		Cimage_header.o dskio.o Cstring.o dtreklogic.o dtrekvec.o \
		Ccrystal.o Cdetector.o Csource.o Cwavelength.o Cstat.o \
		CBeamData.o CBackground.o dtarray.o CScanBitmap.o \
		Cgoniometer.o  dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)


CBackground.o:	CBackground.cc CBackground.h

CBeamData.o:	CBeamData.cc CBeamData.h

dtdisplayupdate.o:	dtdisplayupdate.cc Dtrek.h Cstring.cc

dtdisplayupdate:	dtdisplayupdate.o Cstring.o dtreksys.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtintegrate.o:	dtintegrate.cc Dtrek.h Cintegrate.h Cscan.h Creflnlist.h \
		Crotation.h $(ANL_TIMER_INC)
		$(CC) -DLESSMEM -c -o $@ $(CCFLAGS) $<

dtintegrateh:	dtintegrate.o Cintegrate.o Creflnlist.o Crotation.o Cprofit2D.o \
                 Cspacegroup.o Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                 Cpredict.o Crefine.o Cdetector.o Csource.o Cgoniometer.o \
                 CinterpMesh.o CScanBitmap.o \
                 Ccrystal.o Cimage.o Crefln.o Cimage_header.o Cwavelength.o \
		 Cstring.o Cfind.o dskio.o dtreklogic.o dtrekvec.o \
		 dtreksys.o dtsvd.o dtarray.o dtprofit_main.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

# Pick which dtintegrate[h|l] is the default

dtintegrate:	dtintegrate.o Cintegrate.o Creflnlist.o Crotation.o Cprofit2D.o \
                 Cspacegroup.o Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                 Cpredict.o Crefine.o Cdetector.o Csource.o Cgoniometer.o \
                 Ccrystal.o Cimage.o Crefln.o Cimage_header.o Cwavelength.o \
                 CinterpMesh.o CScanBitmap.o dtarray.o \
		 Cstring.o Cfind.o dskio.o dtreklogic.o dtrekvec.o \
		 dtreksys.o dtsvd.o dtprofit_main.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtprofitv6.o:	dtprofitv6.cc Dtrek.h Cintegrate.h Cscan.h Creflnlist.h \
		Crotation.h

dtprofitv6:	dtprofitv6.o Creflnlist.o Crotation.o Cintegrate.o Cprofit2D.o \
                 Cspacegroup.o Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                 Cpredict.o Crefine.o Cdetector.o Csource.o Cgoniometer.o \
                 Ccrystal.o Cimage.o Crefln.o Cimage_header.o Cwavelength.o \
		 Cstring.o Cfind.o dskio.o dtreklogic.o \
		 dtrekvec.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtprofit_main.o:	dtprofit.cc Dtrek.h Cintegrate.h Cscan.h Creflnlist.h \
		dtprofit.cc Dtrek.h Cintegrate.h Cscan.h Creflnlist.h \
		Crotation.h
		$(CC) -DNO_DTPROFIT_MAIN -c -o $@ $(CCFLAGS) $<

dtprofit.o:	dtprofit.cc Dtrek.h Cintegrate.h Cscan.h Creflnlist.h \
		Crotation.h

dtprofit:	dtprofit.o Creflnlist.o Crotation.o Cintegrate.o Cprofit2D.o \
                 Cspacegroup.o Cscan.o Cnonunf.o BitmapRLE.o Cspatial.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o \
                 Cpredict.o Crefine.o Cdetector.o Csource.o Cgoniometer.o \
                 Ccrystal.o Cimage.o Crefln.o Cimage_header.o Cwavelength.o \
                 CinterpMesh.o CScanBitmap.o dtarray.o \
		 Cstring.o Cfind.o dskio.o dtreklogic.o \
		 dtrekvec.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtimagederiv.o:	dtimagederiv.cc Dtrek.h Cscan.h Cimage.h Cimage_header.h

dtimagederiv:	dtimagederiv.o Cscan.o Cimage.o Cstring.o Crotation.o \
		Cimage_header.o Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o Crefln.o \
		Creflnlist.o Ccrystal.o Cspacegroup.o dtreklogic.o \
		dskio.o dtrekvec.o dtreksys.o dtsvd.o dtarray.o $(MSCOBJS)  
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtadd2.o:	dtadd2.cc Dtrek.h Crotation.h Cimage.h Cimage_header.h

dtadd2:		dtadd2.o Cscan.o Cimage.o Cstring.o Crotation.o \
		Cimage_header.o Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o Crefln.o \
		Creflnlist.o Ccrystal.o Cspacegroup.o Cdetector.o Csource.o \
		Cgoniometer.o Cspatial.o Cwavelength.o dtsvd.o dtarray.o \
		dtreklogic.o dskio.o dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtcompare:	dtcompare.o Ccompare.o Cstring.o Creflnlist.o \
		dtreklogic.o dtreksys.o dskio.o dtrekvec.o \
		Cimage_header.o Crotation.o Cgoniometer.o \
		Cimage.o Crefln.o Cspacegroup.o Ccrystal.o \
		Cstat.o Csource.o Cdetector.o \
		Cscan.o Cspatial.o Cnonunf.o BitmapRLE.o \
		Cwavelength.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o dtsvd.o dtarray.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)
                                                                          
dtcompare.o:    dtcompare.cc Ccompare.h

Ccompare.o:     Ccompare.cc Ccompare.h

dtaverage.o:	dtaverage.cc Dtrek.h Cscan.h Cimage.h Cimage_header.h

dtaverage:	dtaverage.o Cscan.o Cimage.o Cstring.o Crotation.o \
		Cimage_header.o Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o Crefln.o \
		Creflnlist.o Ccrystal.o Cspacegroup.o Cdetector.o Csource.o \
		Cgoniometer.o Cspatial.o Cwavelength.o dtsvd.o dtarray.o \
		dtreklogic.o dskio.o dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtfind.o:	dtfind.cc Dtrek.h Cfind.h Cscan.h Creflnlist.h Crotation.h

dtpredictimage:	dtpredictimage.o Cpredictimage.o Creflnlist.o Cscan.o \
		Cpredict.o Crotation.o  Cnonunf.o BitmapRLE.o Cspatial.o Cspacegroup.o \
                Cgoniometer.o Ccrystal.o Cwavelength.o Csource.o Cdetector.o \
                Cimage.o Crefln.o Cimage_header.o Cstring.o dtsvd.o dtarray.o \
		C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dskio.o dtrekvec.o dtreklogic.o dtreksys.o \
		$(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtpredict:	dtpredict.o Cpredict.o Creflnlist.o Cscan.o Crotation.o \
                Cnonunf.o BitmapRLE.o Cspatial.o Cspacegroup.o Cgoniometer.o\
		Ccrystal.o Cwavelength.o Csource.o Cdetector.o \
                Cimage.o Crefln.o Cimage_header.o Cstring.o dtsvd.o dtarray.o \
		C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dskio.o dtrekvec.o dtreklogic.o dtreksys.o \
		$(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtpredict.o:	dtpredict.cc Dtrek.h Cfind.h Cscan.h Creflnlist.h Crotation.h\
			Ccrystal.h Csource.h Cdetector.h Cgoniometer.h

dtrekvec.o:	dtrekvec.cc dtrekvec.h

dtreklogic.o:	dtreklogic.cc dtreklogic.h

dtreksys.o:	dtreksys.cc dtreksys.h

dtreksys.h:	Cstring.h Dtrek.h

C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o:	C3Ddata.cc C3Ddata.h 

C3Ddata.h:      Dtrek.h Cimage.h

Cstrategy.o:	Cstrategy.cc Cstrategy.h

Cstrategy.h:	Dtrek.h Cimage_header.h

Cpredict.o:	Cpredict.cc Cpredict.h Dtrek.h

Cpredict.h:	Cimage_header.h Cdetector.h Csource.h Ccrystal.h \
		Creflnlist.h dtrekvec.h

Corient.o:	Corient.cc Corient.h Dtrek.h

Corient.h:	Cimage_header.h Csource.h Ccrystal.h Cgoniometer.h \
		Crotation.h dtrekvec.h

Crefine.o:	Crefine.cc Crefine.h Dtrek.h

Crefine.h:	Cimage_header.h Cdetector.h Csource.h Ccrystal.h \
		Crotation.h Creflnlist.h dtrekvec.h

CCellReduce.o:	CCellReduce.cc CCellReduce.h

CCellReduce.h:	Dtrek.h Cimage_header.h Cdetector.h Csource.h Ccrystal.h \
		Crotation.h Creflnlist.h C3Ddata.h Cimage.h dtrekvec.h

CEyeballIndex.o:	CEyeballIndex.cc CEyeballIndex.h Dtrek.h Creflnlist.h Cimage_header.h

Cindex.o:	Cindex.cc Cindex.h Dtrek.h

Cindex.h:	Cimage_header.h Cdetector.h Csource.h Ccrystal.h \
		Crotation.h Creflnlist.h dtrekvec.h

Cfind.o:	Cfind.cc Dtrek.h Cfind.h Cimage.h Cscan.h Creflnlist.h

Cfind.h:	Cimage.h Crotation.h Cspatial.h Cnonunf.h Cscan.h Creflnlist.h

dskio.o:	dskio.cc dskio.h $(ANL_DTREKIO)

raxis.o:	raxis.cc raxis.h

bs.o:		bs.cc bs.h raxis.h

marccd.o:	marccd.cc marccd.h

marip.o:	marip.cc marip.h

winbmp.o:	winbmp.cc winbmp.h

pck.o:		pck.c

adscccd.o:	adscccd.cc adscccd.h

brandeisccd.o:	brandeisccd.cc brandeisccd.h

dip2030.o:	dip2030.cc dip2030.h

medoptics.o:	medoptics.cc medoptics.h

DevErrCode.o:	DevErrCode.cc DevErrCode.h

tCimage_header:	tCimage_header.o Cimage_header.o Cstring.o dskio.o \
		Cimage.o Crefln.o Creflnlist.o Ccrystal.o Cspacegroup.o \
		dtreklogic.o dtrekvec.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tCimage_header.o:	tCimage_header.cc Cimage_header.h

tMSCstring:	tMSCstring.o MSCstring.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) -lm

MSCstring.o:	MSCstring.cc MSCstring.h

Cstring.o:	Cstring.cc Cstring.h

dtrek2raxis:	dtrek2raxis.o Cimage.o Cimage_header.o Cstring.o dskio.o \
		Crefln.o Creflnlist.o Cspacegroup.o dtreklogic.o \
		Cdetector.o Csource.o Cspatial.o Cnonunf.o BitmapRLE.o Cspacegroup.o \
		Cwavelength.o Cgoniometer.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dtsvd.o dtarray.o \
		Crotation.o Cscan.o \
		Ccrystal.o dtrekvec.o  dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtrek2raxis.cc:	Cimage.h

dtraxis:	dtraxis.o Cimage.o Cimage_header.o Cstring.o dskio.o \
		Crefln.o Creflnlist.o Cspacegroup.o dtreklogic.o \
		Ccrystal.o dtrekvec.o  dtreksys.o Cspatial.o Csource.o \
		Cdetector.o Cgoniometer.o Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dtsvd.o dtarray.o \
		Crotation.o Cscan.o \
		Cwavelength.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtraxis.cc:	Cimage.h

dtreorder:	dtreorder.o Cimage.o Cimage_header.o Cstring.o Crefln.o \
		Ccrystal.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtnonunfedit:	dtnonunfedit.o Cimage.o Cimage_header.o Cstring.o Crefln.o \
		Cdetector.o Cspatial.o Cnonunf.o BitmapRLE.o Csource.o Cwavelength.o \
		Cgoniometer.o Cstring.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dtsvd.o dtarray.o \
		Ccrystal.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		dtrekvec.o dtreksys.o Cnonunf.o BitmapRLE.o Crotation.o Cscan.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtflip:		dtflip.o Cimage.o Cimage_header.o Cstring.o Crefln.o \
		Ccrystal.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtsmooth:	dtsmooth.o Cimage.o Cimage_header.o Cstring.o Crefln.o dtsvd.o dtarray.o \
		Ccrystal.o Creflnlist.o Cspacegroup.o dskio.o dtreklogic.o \
		Cgoniometer.o Csource.o Cwavelength.o Cdetector.o Cspatial.o \
		Crotation.o Cscan.o \
		Cnonunf.o BitmapRLE.o C3Ddata.o CScanBitmap.o Cprofit2D.o Cgenetic.o CEmpericalProfit.o CinterpMesh.o Cstat.o dtrekvec.o dtreksys.o $(MSCOBJS)
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

dtsvd.o:	dtsvd.h dtsvd.cc

dtarray.o:	dtarray.h dtarray.cc

install:	$(ALL)
		if [ ! -d $(INSTALLDIR) ];then mkdir $(INSTALLDIR);fi;
		for X in $(ALL); do \
			cp $$X $(INSTALLDIR); done

clean:
		$(RM) *.o *$(OSTYPE).a $(ALL)


cleanbin:
		$(RM) $(ALL) Cimage_header.o Cimage.o Cstring.o Creflnlist.o dt*.o

cleancpp:
		 ls -la | awk '/.*\.cpp/{ print thad $$9 ; }' | \
		awk -F. '{ print "rm " $$1".cc" ;}' | sh

tTests:		tTests.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tTests.o:	tTests.cc

tio:		tio.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tio.o:		tio.cc

tss:		tss.o
		$(CC) -o $@ $^  $(CCFLAGS) $(GLIBDIR) $(LIBS)

tss.o:		tss.cc

.cc.o:
		$(COMPILE.cc) $< 

.c.o:
		$(COMPILE.cc) $< 
8mm:
		tar cf - /rxdat1/people/jwp/DTREK* \
			/user5/jwp/DTREK* \
		 | rsh guest@tapes dd of=/dev/rst0 obs=1024

8mm.tarlis:
		rsh guest@tapes dd if=/dev/rst0 ibs=1024 | tar tvf - >& $@

