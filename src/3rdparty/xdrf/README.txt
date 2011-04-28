libxdrf version 1.1

This directory contains source of the portable data compression 
library xdrf, which was developed for EUROPORT.
It also contains two test programs (one in written in FORTRAN and
one written in C) and a large file with test data.

To create the library, look into the subdirectory conf, and select
the architecture (ARCH) that best matches your system. Edit the makefile
and use your system definition.
Note that you have to set ARCH, LIBS, HASRANLIB and perhaps some
compiler specific settings.

Type;

make

This will create libxdrf.a and the two test programs (both programs
are linked with this new libxdrf.a library).

To test the program type;
ctest 	(FORTRAN version), or 
ftest   (C version)

You can compare the input file 'test.gmx' with the output 'test.out',
which should be the same (numerically at least).

You can also compare the filesize of 'test.gmx' and the newly-created
compressed data file 'test.xdr'. 

Read the included manual page, or look into the source files.
'Intro.txt' describes in general terms what xdrf is for.

- frans

(c) 1995 frans van hoesel

hoesel@chem.rug.nl

