# Master Makefile for d*TREK 
#
# Copyright (C) 2017 Rigaku Americas Corporation
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

# Edit the next lines for your site
#
# Directory where OpenMotif is installed.
MOTIFHOME  := /Applications/OpenMotif21
# Directory where CBF library is located.
CBFHOME    := $(CURDIR)/CBF
# Directory where PTypes library is located.
PTYPESHOME := $(CURDIR)/ptypes

# Lines below should not need to be modified

$(info )
$(info d*TREK master Makefile)

ifeq (,$(PTYPESHOME))
  $(error PTYPESHOME not found - edit Makefile with correct directory path)
endif
ifeq (,$(CBFHOME))
  $(error CBFHOME not found - edit Makefile with correct directory path)
endif

SHELL := /bin/sh
OS    := $(shell [ "Darwin" == "`uname`" ] && echo osx || echo linux)
ARCH  := $(shell uname -m)

$(info )
$(info Operating system is $(OS))
$(info Architecture is $(ARCH))

OSTYPE_EXT.osx   = darwin11
OSTYPE_EXT.linux = linux
GMAKE.osx        = /usr/bin/make
GMAKE.linux      = /usr/bin/gmake

OSTYPE     = _$(OSTYPE_EXT.$(OS))
GMAKE      = $(GMAKE.$(OS))
MFILE      = Makefile.$(OS)
INSTALLDIR = $(CURDIR)/bin$(OSTYPE)

CC  = g++
CXX = g++

$(info )
$(info Motif directory is $(MOTIFHOME))
$(info CBF library directory is $(CBFHOME))
$(info PTypes library directory is $(PTYPESHOME))
$(info Installation directory is $(INSTALLDIR))

export OSTYPE
export CC
export CXX
export MOTIFHOME
export CBFHOME
export PTYPESHOME
export INSTALLDIR

$(info )
all:
	cd DTTREK;    $(GMAKE) -f $(MFILE) lib
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) lib
	cd DTDEV;     $(GMAKE) -f $(MFILE) lib
	cd DTDISPLAY; $(GMAKE) -f $(MFILE)
	cd DTTREK;    $(GMAKE) -f $(MFILE)
	cd DTPROCESS; $(GMAKE) -f $(MFILE)

clean:
	cd DTTREK;    $(GMAKE) -f $(MFILE) clean
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) clean
	cd DTDEV;     $(GMAKE) -f $(MFILE) clean
	cd DTPROCESS; $(GMAKE) -f $(MFILE) clean

cleanbin:
	cd DTTREK;    $(GMAKE) -f $(MFILE) cleanbin
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) cleanbin
	cd DTPROCESS; $(GMAKE) -f $(MFILE) cleanbin

install:
	cd DTTREK;    $(GMAKE) -f $(MFILE) install
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) install
	cd DTPROCESS; $(GMAKE) -f $(MFILE) install

