# Introduction #

The d\*TREK command line programs are built from the command line using make.  Depending on the operating system, there are specific makefiles that should be used; generally these are named Makefile.$os, where $os is the name of the operating system (Makefile.linux for Linux, etc).


# Linux #

Building the command line programs on Linux is a two step process:
  * In the DTDISPLAY subdirectory, run the command
```
       make -f Makefile.linux lib
```
  * In the DTTREK subdirectory, run the command
```
       make -f Makefile.linux
```