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
      SUBROUTINE WRITE_LAPS(REFTIME,VALTIME,DIR,EXT,IMAX,JMAX,
     1   KMAX,KDIM,VAR,LVL,LVL_COORD,UNITS,COMMENT,DATA,
     1   ISTATUS)

C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            writelapsdata
C            cvt_str_data
C
C      The writelapsdata subroutine reads the following FORTRAN
C      subroutines from the readlapsdata.for file:
C
C            cvt_fname_data
C
C      The writelapsdata subroutine reads the following C subroutines
C      from the readwritelaps.c file:
C            make_c_fname
C            write_cdf_file
C            cre_lw3, cre_lh1, cre_lh2, cre_lh3, cre_lh4, cre_lq3,
C            cre_lsx, cre_lwm, cre_lt1, cre_lhe, cre_liw, cre_lmt,
C            cre_lmr, cre_lf1, cre_l1s, cre_lps, cre_lrp, cre_lba,
C            cre_lc3, cre_lwc, cre_lil, cre_lcb, cre_lct, cre_lcv,
C            cre_lmd, cre_lco, cre_lty, cre_lcp
C            write_header_cdf
C            cdf_update_laps
C            cdf_get_index
C            cdf_get_coord
C            cdf_dim_size
C            cdf_write_grid
C            cdf_update_laps_inv
C
C**********************************************************************
C
C      Subroutine WRITE_LAPS_DATA
C
C      Author:    John Snook
C      Modified:  To write netCDF data files      1/93 Linda Wharton
C                 To remove BYTE arrays           4/94 Linda Wharton
C                 To remove downcase of filename  9/95 Linda Wharton
C
C      Writes data in arrays DATA and COMMENT to the netCDF file name
C      specified by I4TIME, DIR and EXT.  The data in VAR, LVL, LVL_COORD,
C      IMAX, JMAX, KMAX, KDIM and UNITS are stored into the netCDF file
C      when it is created.  ISTATUS is returned.
C
C**********************************************************************
C
        IMPLICIT        NONE
C
        INTEGER*4        ML
        PARAMETER       (ML=10000)

        INTEGER*2       FN_LENGTH
C
        INTEGER*4       REFTIME,         !model runtime of data
     1          VALTIME,         !valid time of data
     1          IMAX,JMAX,KMAX, !# cols, # rows, # fields
     1          KDIM,           !K dimension of DATA array
     1          LVL(*),         !Level of each field (4 digit max)
     1          FLAG,
     1          ERROR(2),
     1          I,ISTATUS
C
        REAL*4          DATA(IMAX,JMAX,KDIM)    !Raw data to be written
C
        CHARACTER*(*)    DIR             !Directory to be written to
        CHARACTER*31    EXT             !File name ext (up to 31 chars)
        CHARACTER*31    EXT_IN          !INPUT File name ext (up to 31 chars)
        CHARACTER*3     VAR(*)          !3 letter ID of each field
        CHARACTER*3     VAR_IN(300)     !INPUT 3 letter ID of requested fields
        CHARACTER*4     LVL_COORD(*)    !Vertical coordinate for each field
        CHARACTER*10    UNITS(*)        !units of each field
        CHARACTER*125   COMMENT(*)      !Comments for each field
        CHARACTER*9     GTIME
        CHARACTER*4     fcst_hh_mm
        CHARACTER*91    FILE_NAME
        CHARACTER*11    LAPS_DOM_FILE   !Name of domain file e.g. NEST7GRID
        CHARACTER*24    ASCTIME
C
        COMMON          /PRT/FLAG
C
C-------------------------------------------------------------------------------
C
        ERROR(1)=1
        ERROR(2)=0
        ext_in = ext
C
C ****  Various checks on input data.
C
        IF (KMAX .GT. KDIM) THEN
                IF (FLAG .NE. 1)
     1   write (6,*) 'Illegal K dimension in DATA array...write aborted.
     1'
                ISTATUS=ERROR(2)
                RETURN
        ENDIF
C
C ****  Create ascii time variables.
C
        CALL CV_I4TIM_ASC_LP(VALTIME,ASCTIME,ISTATUS)
C
C ****  Specify LAPS_DOM_FILE
C
        LAPS_DOM_FILE = 'nest7grid'
C
C ****  Specify file name
C
        CALL MAKE_FNAM_LP(REFTIME,GTIME,ISTATUS)
        IF (ISTATUS .ne. 1) THEN
!               CALL LOG_ERROR_GG(' ',ISTATUS,
!       1         'Error converting i4time to file name...save aborted.',0,0)
                ISTATUS=ERROR(2)
                RETURN
        ENDIF

        call make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)

        do i=1,kdim
          call upcase(var(i),var_in(i))
        enddo

        call cvt_fname_data(dir,gtime,fcst_hh_mm,ext,
     1                    file_name,fn_length,istatus)
        call upcase(ext,ext_in)

C
C **** write out netCDF file
C
      CALL WRITE_CDF_FILE(FILE_NAME,FN_LENGTH,LAPS_DOM_FILE,
     1                   ASCTIME,EXT_IN,VAR_IN,LVL_COORD,UNITS,
     1                   COMMENT,IMAX,JMAX,KMAX,KDIM,LVL,DATA,ISTATUS)
C
      if (ISTATUS .gt. 0) goto 980
      IF (ISTATUS .eq. -2) GOTO 940
      IF (ISTATUS .eq. -3) GOTO 950
      IF (ISTATUS .eq. -4) GOTO 960
      IF (ISTATUS .eq. -5) GOTO 970
C
C ****  Return normally.
C
        ISTATUS=ERROR(1)
999     RETURN
C
C ****  Error trapping.
C
940     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file to be written to...write abort
     1ed.'
        ISTATUS=ERROR(2)
        GOTO 999
C
950     IF (FLAG .NE. 1)
     1    write (6,*) 'Error in imax,jmax,kmax, or kdim ...write aborted
     1.'
        ISTATUS=ERROR(2)
        GOTO 999
C
960     IF (FLAG .NE. 1)
     1    write (6,*) 'Error writing data to file...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
970     IF (FLAG .NE. 1)
     1    write (6,*) 'Error writing header info into file...write abort
     1ed.'
        ISTATUS=ERROR(2)
        GOTO 999
C
980     IF (FLAG .NE. 1)
     1    write (6,*) 'Some grids not written....could not convert LAPS
     1variables.'
        ISTATUS=ERROR(2)
        GOTO 999
C
        END

C##########################################################################

