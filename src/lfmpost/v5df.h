! v5df.h


! Include file for using v5d functions from FORTRAN programs


! Function prototypes.  See the README file for details.  These are
! the functions you'll want to use for writing v5d file converters.

      integer v5dcreate

      integer v5dcreatesimple

      integer v5dwrite

      integer v5dmcfile

      integer v5dclose


! 5-D grid limits, must match those in v5d.h!!!
      integer MAXVARS, MAXTIMES, MAXROWS, MAXCOLUMNS, MAXLEVELS

      parameter (MAXVARS=400)
      parameter (MAXTIMES=75)
      parameter (MAXROWS=500)
      parameter (MAXCOLUMNS=500)
      parameter (MAXLEVELS=100)

! Missing values
      real MISSING
      integer IMISSING

      parameter (MISSING=1.0E35)
      parameter (IMISSING=-987654)

