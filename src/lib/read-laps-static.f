cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
C Read time-invariant LAPS data.  i4time is set to zero which results
C in a "file time" of Jan 1, 1960 or 1970.
C
        SUBROUTINE READ_LAPS_STATIC (DIR, EXT, IDIM, JDIM, KMAX, KDIM,
     1                          VAR, UNITS, COMMENT, GRID, STATUS)
C
        IMPLICIT NONE
C
        CHARACTER*(*) DIR, EXT
        INTEGER IDIM, JDIM, KMAX, KDIM
        CHARACTER*(*) VAR(1), UNITS(1), COMMENT(1)
        REAL GRID(IDIM,JDIM,KDIM)
        INTEGER STATUS
C
        INTEGER NVARSMAX, I
        PARAMETER (NVARSMAX=20)         ! may need to increase this someday
        CHARACTER*4 LVL_COORD(NVARSMAX)
        INTEGER LVL(NVARSMAX)

        write(6,*)
     1 ' WARNING: This routine reads static data in the old format.'
        write(6,*)
     1 ' The format is non-netCDF and is not portable.'
        write(6,*)
     1 ' It is recommended to call rd_laps_static instead to read'
        write(6,*)
     1 ' the new, portable netCDF files.'

        DO I=1,NVARSMAX
                LVL(I) = 0
        ENDDO

        CALL READ_LAPS_DATA (0, DIR, EXT, IDIM, JDIM, KMAX, KDIM,
     1  VAR, LVL, LVL_COORD, UNITS, COMMENT, GRID, STATUS)

        RETURN
        END
