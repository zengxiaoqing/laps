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
      SUBROUTINE WRITE_LAPS_DATA(I4TIME,DIR,EXT,IMAX,JMAX,
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
        INTEGER*4       I4TIME,         !I4time of data
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
        CALL CV_I4TIM_ASC_LP(I4TIME,ASCTIME,ISTATUS)
C
C ****  Specify LAPS_DOM_FILE
C
        LAPS_DOM_FILE = 'nest7grid'
C
C ****  Specify file name
C
        CALL MAKE_FNAM_LP(I4TIME,GTIME,ISTATUS)
        IF (ISTATUS .ne. 1) THEN
!               CALL LOG_ERROR_GG(' ',ISTATUS,
!       1         'Error converting i4time to file name...save aborted.',0,0)
                ISTATUS=ERROR(2)
                RETURN
        ENDIF
        fcst_hh_mm = '0000'

        do i=1,kdim
          call upcase(var(i),var_in(i))
        enddo

        call cvt_fname_data(dir,gtime,fcst_hh_mm,ext,file_name,
     1                      fn_length,istatus)

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
        SUBROUTINE WRITE_OLD_LAPS(I4TIME,DIR,EXT,IMAX,JMAX,KMAX,KDIM,
     1                     VAR,LVL,LVL_COORD,UNITS,COMMENT,
     1                     DATA,ISTATUS)
C
        IMPLICIT        NONE
C
        INTEGER*4        ML
        PARAMETER       (ML=10000)
C
        INTEGER*4       I4TIME,         !I4time of data
     1          IMAX,JMAX,KMAX, !# cols, # rows, # fields
     1          KDIM,           !K dimension of DATA array
     1          LVL(KDIM),      !Level of each field (4 digit max)
     1          FLAG,
     1          REC,
     1          RCDL,
     1          START_REC,
     1          ERROR(2),
     1          INDX,
     1          I,J,K,
     1          ISTATUS

        INTEGER*2       FN_LENGTH
C
        REAL*4          DATA(IMAX,JMAX,KDIM)    !Raw data to be written
C
        CHARACTER*50    DIR             !Directory to be written to
        CHARACTER*31    EXT             !File name ext (up to 31 chars)
        CHARACTER*3     VAR(KDIM)       !3 letter ID of each field
        CHARACTER*4     LVL_COORD(KDIM) !Vertical coordinate for each field
        CHARACTER*10    UNITS(KDIM)     !units of each field
        CHARACTER*125   COMMENT(KDIM)   !Comments for each field
        CHARACTER*4     CIMAX,CJMAX,CKMAX
        CHARACTER*4     CLVL(ML)
        CHARACTER*4     CSTART_REC
        CHARACTER*9     GTIME
        CHARACTER*4     fcst_hh_mm
        CHARACTER*91    FILE_NAME
        CHARACTER*24    ASCTIME
        CHARACTER*4     MARK
        CHARACTER*1     LETTER(52)
C

        DATA            LETTER/'A','B','C','D','E','F','G','H','I','J',
     1                 'K','L','M','N','O','P','Q','R','S','T',
     1                 'U','V','W','X','Y','Z',
     1                 'a','b','c','d','e','f','g','h','i','j',
     1                 'k','l','m','n','o','p','q','r','s','t',
     1                 'u','v','w','x','y','z'/
C
        COMMON          /PRT/FLAG
C
C-------------------------------------------------------------------------------
C
        ERROR(1)=1
        ERROR(2)=0
C
C ****  Various checks on input data.
C
        IF (IMAX .LT. 1 .OR. IMAX .GT. 9999 .OR.
     1    JMAX .LT. 1 .OR. JMAX .GT. 9999 .OR.
     1    KMAX .LT. 1 .OR. KMAX .GT. 9999) THEN
           IF (FLAG .NE. 1)
     1    write (6,*) 'Illegal value contained in max array dimensions..
     1.',
     1         'Write aborted.'
           ISTATUS=ERROR(2)
           RETURN
        ENDIF
C
        DO K=1,KMAX
                IF (LVL(K) .LT. -999 .OR. LVL(K) .GT. 9999) THEN
                        IF (FLAG .NE. 1)
     1  write (6,*) 'Illegal value contained in level array...Write',
     1                 ' aborted.'
                        ISTATUS=ERROR(2)
                        RETURN
                ENDIF
        ENDDO
C
        IF (KMAX .GT. KDIM) THEN
                IF (FLAG .NE. 1)
     1   write (6,*) 'Illegal K dimension in DATA array...write aborted.
     1'
                ISTATUS=ERROR(2)
                RETURN
        ENDIF
C
C ****  Create character variables of input numerical variables.
C
        write(CIMAX,900)  IMAX
        write(CJMAX,900)  JMAX
        write(CKMAX,900)  KMAX
900     FORMAT(I4)
C
C ****  Create ascii time variables.
C
        CALL CV_I4TIM_ASC_LP(I4TIME,ASCTIME,ISTATUS)
C
C ****  Specify file name and open file.
C
        CALL MAKE_FNAM_LP(I4TIME,GTIME,ISTATUS)
        IF (ISTATUS .ne. 1) THEN
!               CALL LOG_ERROR_GG(' ',ISTATUS,
!       1         'Error converting i4time to file name...save aborted.',0,0)
                ISTATUS=ERROR(2)
                RETURN
        ENDIF

        fcst_hh_mm = '0000'
        call cvt_fname_data(dir,gtime,fcst_hh_mm,ext,file_name,
     1                      fn_length,istatus)
        RCDL=IMAX+1
        IF (RCDL .LT. 40) RCDL=40
        OPEN(1,FILE=FILE_NAME,STATUS='NEW',
     1     ACCESS='DIRECT',RECL=RCDL,ERR=940)
C
C ****  Write header record number 1.  Contains:
C          1)  Record length (IMAX),
C          2)  Number of lines (JMAX),
C          3)  Number of fields (KMAX),
C          4)  Ascii time,
C          5)  Version number of write routine.
C
        REC=1
        WRITE(1,REC=REC) '*   ',CIMAX,CJMAX,CKMAX,ASCTIME(1:17),'  V1'
C
C ****  Write one header record for each field containing:
C          1)  Variable (VAR),
C          2)  Level (LVL),
C          3)  Vertical coordinate (LVL_COORD)
C          4)  Units,
C          5)  Starting record number of data,
C          6)  Comments.
C
        START_REC=2+KMAX
        DO K=1,KMAX
                REC=REC+1
                IF (START_REC .GE. 62000) THEN
                print *,'Number of data records exceeds 62,000...Write a
     1borted.'
                        ISTATUS=ERROR(2)
                        RETURN
                ENDIF
                IF (START_REC .LT. 10000) THEN
                        write(CSTART_REC,900)  START_REC
                ELSE
                        INDX=START_REC/1000-9
                        write(CSTART_REC(2:4),901)  MOD(START_REC,1000)
                        IF (CSTART_REC(2:2) .EQ. ' ') CSTART_REC(2:2)='0
     1'
                        IF (CSTART_REC(3:3) .EQ. ' ') CSTART_REC(3:3)='0
     1'
                        CSTART_REC(1:1)=LETTER(INDX)
                ENDIF
                write(CLVL(K),900)  LVL(K)
                CALL UPCASE(VAR(K),VAR(K))
                CALL UPCASE(LVL_COORD(K),LVL_COORD(K))
                CALL UPCASE(UNITS(K),UNITS(K))
                WRITE(1,REC=REC) '*',VAR(K),CLVL(K),LVL_COORD(K),
     1                               UNITS(K),CSTART_REC,COMMENT(K)
                START_REC=START_REC+JMAX
        ENDDO
901     FORMAT(I3)
C
C ****  Write data.
C
        DO K=1,KMAX
           IF (LVL(K) .LT. 100) THEN
                MARK=VAR(K)(1:2)//CLVL(K)(3:4)
           ELSEIF (LVL(K) .LT. 1000) THEN
                MARK=VAR(K)(1:2)//CLVL(K)(2:3)
           ELSE
                MARK=VAR(K)(1:2)//CLVL(K)(1:2)
           ENDIF
           DO J=1,JMAX
                REC=REC+1
                WRITE(1,REC=REC) MARK,(DATA(I,J,K),I=1,IMAX)
           ENDDO
        ENDDO
C
C ****  Return normally.
C
        IF (FLAG .NE. 1) write (6,*) REC,' records written.'
        ISTATUS=ERROR(1)
998     CLOSE(1,ERR=999)
999     RETURN
C
C ****  Error trapping.
C
940     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file to be written to...write abort
     1ed.'
        ISTATUS=ERROR(2)
        GOTO 998
C
        END
C##########################################################################

