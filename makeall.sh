#!/bin/bash

[ "Darwin" = "`uname`" ] && MAKE=make || MAKE=gmake
[ "Darwin" = "`uname`" ] && BINDIR=bin_darwin11 || BINDIR=bin_linux
[ "Darwin" = "`uname`" ] && OS=osx || OS=linux

homeDir=`dirname $0`
[ ! "\." = "${homeDir}" ] && cd $homeDir

echo
echo Running ${MAKE} in $PWD
echo using Makefile $OS variants
echo placing results in $BINDIR
echo

${MAKE} -f Makefile clean
level=$? && [ $level -ne 0 ] && exit $level
${MAKE} -f Makefile
level=$? && [ $level -ne 0 ] && exit $level
${MAKE} -f Makefile install
level=$? && [ $level -ne 0 ] && exit $level


exit 0

