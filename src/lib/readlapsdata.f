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
      SUBROUTINE READ_LAPS_DATA(I4TIME,DIR,EXT,IIMAX,JJMAX,KKMAX,KDIM,
     1                    VAR_REQ,LVL_REQ,LVL_COORD_REQ,
     1                    UNITS_REQ,COMMENT_REQ,DATA,ISTATUS)
C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            readlapsdata
C            get_packed_data
C            cvt_fname_data
C            cvt_var_data
C            cvt_hdr_char
C            cvt_byte_char
C
C      The readlapsdata subroutine reads the following C subroutines
C      from the readwritelaps.c file:
C            make_c_fname
C            read_cdf_file
C            get_cdf_var
C            open_cdf
C            cdf_retrieve_header
C            itoa
C            cdf_retrieve_laps_grid
C            cdf_get_index
C            cdf_get_coord
C            cdf_dim_size
C            cdf_chk_laps_inv
C            fill_empty_grids
C
C**********************************************************************
C
C      Subroutine READ_LAPS_DATA
C
C      Author:    John Snook
C      Modified:  To accept netCDF data files  1/93 Linda Wharton
C
C      Reads data requested by arrays VAR_REQ and LVL_REQ for the
C      I4TIME, DIR and EXT specified.  Returns LVL_COORD-REQ,
C      UNITS_REQ, COMMENT_REQ, DATA AND ISTATUS
C
C**********************************************************************
C
      IMPLICIT  NONE
C
      INTEGER*4       MF
      PARAMETER (MF=300)                !Max fields allowed

      INTEGER*2 FN_LENGTH

      INTEGER*4 I4TIME,                 !INPUT I4time of data
     1          IIMAX,JJMAX,KKMAX,      !INPUT # cols, # rows, # fields
     1          KDIM,                   !INPUT K dimension of DATA array
     1          LVL_REQ(KDIM),          !INPUT Requested levels
     1          FLAG,                   !Print flag (1 = off)
     1          no_laps_diag,           !if = 0, print diagnostic output
     1          REC,
     1          RCDL,
     1          IMAX,JMAX,KMAX,
     1          LVL(MF),
     1          START_REC(MF),
     1          I,J,K,KK,
     1          ERROR(3),
     1          ISTATUS                     !OUTPUT
C
      Integer*4 process_var(300),num_variables
C
      REAL*4    DATA(IIMAX,JJMAX,KDIM),     !OUTPUT data
     1          MSG_FLAG
C
      CHARACTER*(*)      DIR                 !INPUT Directory to read data from
      CHARACTER*31      EXT                 !INPUT File name ext (up to 31 chars)
      CHARACTER*31      EXT_IN              !INPUT File name ext (up to 31 chars)
      CHARACTER*3       VAR_REQ(KDIM)       !INPUT 3 letter ID of requested fields
      CHARACTER*3       VAR_IN(MF)          !INPUT 3 letter ID of requested fields
      CHARACTER*4       LVL_COORD_REQ(KDIM) !OUTPUT Vertical coordinate of fields
      CHARACTER*10      UNITS_REQ(KDIM)     !OUTPUT Units of requested fields
      CHARACTER*125     COMMENT_REQ(KDIM)   !OUTPUT Comments for requested fields
      CHARACTER*9       GTIME
      CHARACTER*4       fcst_hh_mm
      CHARACTER*91      FILE_NAME
      CHARACTER*4       MARK
      CHARACTER*4       CIMAX,CJMAX,CKMAX
      CHARACTER*11      LAPS_DOM_FILE   !Name of domain file e.g. NEST7GRID
      character*12      LDF             !LAPS_DOM_FILE + 1 for tfr to C
      CHARACTER*24      ASCTIME
      character*25      asct            !ASCTIME + 1 for tfr to C
      CHARACTER*4       VERSION
      CHARACTER*131     MODEL           !Meteorological model in file
      CHARACTER*131     ORIGIN          !Location where file was created
      CHARACTER*1       HMARK
      CHARACTER*3       VAR(MF)
      CHARACTER*4       CLVL
      CHARACTER*5       LVL_COORD(MF)
      CHARACTER*11      UNITS(MF)
      CHARACTER*4       CSTART_REC
      CHARACTER*126     COMMENT(MF)
C
      LOGICAL         l_packed_data,l_some_missing
C
      DATA              MSG_FLAG/9.E30/
C
      COMMON            /PRT/FLAG
      COMMON            /laps_diag/no_laps_diag
      COMMON            /LAPS_P_VAR/process_var
C
C-------------------------------------------------------------------------------
C
      ERROR(1)=1
      ERROR(2)=0
      ERROR(3)=-2
      l_some_missing = .false.
      ext_in = ext
C
C ****  Create file name.
C
      CALL MAKE_FNAM_LP(I4TIME,GTIME,ISTATUS)
      IF (ISTATUS .ne. 1) THEN
!               CALL LOG_ERROR_GG(' ',ISTATUS,
!       1         'Error converting i4time to file name...read aborted.',0,0)
                ISTATUS=ERROR(2)
                RETURN
      ENDIF
      fcst_hh_mm = '0000'

      do i=1,kdim
        call upcase(var_req(i),var_in(i))
        process_var(i) = 0
      enddo

      call cvt_fname_data(dir,gtime,fcst_hh_mm,ext,file_name,
     1                    fn_length,istatus)

      call upcase(ext,ext_in)

      call read_cdf_file(iimax,jjmax,kkmax,kdim,imax,jmax,kmax,
     1                   num_variables,var_in,lvl_req,file_name,
     1                   fn_length,ext_in,lvl_coord,units,
     1                   comment,ldf,asct,version,model,origin,data,
     1                   process_var,no_laps_diag,istatus)

      if (istatus .ge. 0) then   !return from read with no errors
                                 !convert byte data to characters

         do i = 1, kdim
           lvl_coord_req(i) = lvl_coord(i)
           units_req(i) = units(i)
           comment_req(i) = comment(i)
         enddo
         laps_dom_file = ldf
         asctime = asct

         if (istatus .gt. 0) l_some_missing = .true.

         if(l_some_missing) then
            ISTATUS=ERROR(3)
         else
            ISTATUS=ERROR(1)
         endif
         goto 999
      endif


      if (istatus .eq. -3)  goto 980  !error in dimensioning arrays
      if (istatus .eq. -2)  goto 970  !error retrieving data
      if (istatus .eq. -1) then

C  error opening file as netCDF, try to
C  read file in old LAPS format

      if (no_laps_diag .eq. 0)
     1         write(6,*) 'Cannot open netCDF file - ',
     1                    'Reading file in old LAPS format'


C      IF (ext(1:3) .eq. 'LC3' .or. ext(1:3) .eq. 'SC3' .or.
C     1    ext(1:3) .eq. 'LPS' .or. ext(1:3) .eq. 'SPS' .or.
C     1    ext(1:3) .eq. 'LRP' .or. ext(1:3) .eq. 'SRP' .or.
C     1    ext(1:3) .eq. 'LTY' .or. ext(1:3) .eq. 'STY' .or.
C     1    ext(1:3) .eq. 'LH3' .or. ext(1:3) .eq. 'SH3' .or.
C     1    ext(1:3) .eq. 'LMD' .or. ext(1:3) .eq. 'SMD' )then
        l_packed_data = .false.
C      else
C        l_packed_data = .false.
C      endif

C
C ****  Open file and read first header record.
C
      OPEN(1,FILE=FILE_NAME,STATUS='OLD',ERR=940)

      if(l_packed_data) then
        CLOSE(1)
        goto1000
      endif

        READ(1,890,ERR=895) MARK,CIMAX,CJMAX,CKMAX,ASCTIME,VERSION
890     FORMAT(A4,3A4,A17,A4)
        CLOSE(1)
        IF (MARK(1:1) .NE. '*') THEN
          goto 895
        ENDIF
C
C ****  Decode charater variables.
C
        read(CIMAX,900)  IMAX
        read(CJMAX,900)  JMAX
        read(CKMAX,900)  KMAX
900     FORMAT(I4)
C
        IF (IMAX .GT. IIMAX .OR. JMAX .GT. JJMAX .OR. KKMAX .GT. KDIM) T
     1HEN
                IF (FLAG .NE. 1)
     1          write (6,*) 'Data array dimensioned too small...read abo
     1rted.'
                ISTATUS=ERROR(2)
                GOTO 998
        ENDIF
C
C ****  Open file for direct access.
C
        RCDL=IMAX+1
        IF (RCDL .LT. 40) RCDL=40
        OPEN(1,FILE=FILE_NAME,STATUS='OLD',
     1     ACCESS='DIRECT',
     1     RECL=RCDL,ERR=950)
C
C ****  Read header files for each field.
C
        REC=1
        DO K=1,KMAX
                REC=REC+1
                READ(1,REC=REC,ERR=905) HMARK,VAR(K),CLVL,LVL_COORD(K),U
     1NITS(K),
     1              CSTART_REC,COMMENT(K)
                IF (HMARK .NE. '*') THEN
                  goto 905
                ENDIF
                read(CLVL,900)  LVL(K)
                IF (ICHAR(CSTART_REC(1:1)) .GE. 65 .AND.
     1      ICHAR(CSTART_REC(1:1)) .LE. 90) THEN
                        read(CSTART_REC(2:4),906)  START_REC(K)
                        START_REC(K)=START_REC(K)+10000
     1                      +FLOAT(ICHAR(CSTART_REC(1:1))-65)*1000.
                ELSEIF (ICHAR(CSTART_REC(1:1)) .GE.  97 .AND.
     1          ICHAR(CSTART_REC(1:1)) .LE. 122) THEN
                        read(CSTART_REC(2:4),906)  START_REC(K)
                        START_REC(K)=START_REC(K)+36000
     1                      +FLOAT(ICHAR(CSTART_REC(1:1))-97)*1000.
                ELSE
                        read(CSTART_REC,900)  START_REC(K)
                ENDIF
        ENDDO
906     FORMAT(I3)
C
C ****  Fill data array for requested data.
C
1000    DO 100 KK=1,KKMAX

           if(l_packed_data)then
             CALL UPCASE(VAR_REQ(KK),VAR_REQ(KK))
C             call get_packed_data(i4time,ext,VAR_REQ(kk),LVL_REQ(KK)
C     1             ,DATA(1,1,KK),IIMAX,JJMAX,gtime,file_name,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' Error in get_packed_data'
                 return
             endif

           else
             CALL UPCASE(VAR_REQ(KK),VAR_REQ(KK))
910          FORMAT(1X,A3,I4,' field not found.')
             DO K=1,KMAX
                IF (VAR_REQ(KK) .EQ. VAR(K) .AND. LVL_REQ(KK) .EQ. LVL(K
     1)) THEN
                   REC=START_REC(K)
                   LVL_COORD_REQ(KK)=LVL_COORD(K)
                   UNITS_REQ(KK)=UNITS(K)
                   COMMENT_REQ(KK)=COMMENT(K)
                   DO J=1,JMAX
                        READ(1,REC=REC,ERR=960) MARK,(DATA(I,J,KK),I=1,I
     1MAX)
                        REC=REC+1
                   ENDDO
                   GOTO 100
                ENDIF
             ENDDO
             write(6,910)VAR_REQ(KK),LVL_REQ(KK)
             DO J=1,JMAX
             DO I=1,IMAX
                DATA(I,J,KK)=MSG_FLAG
             ENDDO
             ENDDO
             IF (KKMAX .EQ. 1) THEN
                l_some_missing = .true.
                GOTO 998
             ELSE
                l_some_missing = .true.
             ENDIF

          endif ! Packed Data

100     CONTINUE

      endif !if file fails to open as netcdf
C
C ****  Return normally.
C
        ISTATUS=ERROR(1)
998     if(l_some_missing)ISTATUS=ERROR(3)
        CLOSE(1,ERR=999)
999     RETURN
C
C ****  Error trapping.
C
895     IF (FLAG .NE. 1)
     1    write (6,*) 'Error reading main header...read aborted.'
        ISTATUS=ERROR(2)
        GOTO 998
C
905     IF (FLAG .NE. 1)
     1    write (6,*) 'Error reading field header...read aborted.'
        ISTATUS=ERROR(2)
        GOTO 998
C
940     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file for given i4time...read aborte
     1d.'
        ISTATUS=ERROR(2)
        GOTO 998
C
950     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file as direct access...read aborte
     1d.'
        ISTATUS=ERROR(2)
        GOTO 998
C
960     IF (FLAG .NE. 1)
     1    write (6,*) 'Error finding record number...read aborted.'
        ISTATUS=ERROR(2)
        GOTO 998
C
970     IF (FLAG .NE. 1)
     1    write (6,*) 'Error retrieving data...read aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
980     IF (FLAG .NE. 1)
     1    write (6,*) 'Error in array dimensioning...read aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
        END

C########################################################################
      subroutine cvt_fname_data(dir,rtime,fhh,ext,file_name,
     1                          fn_length,istatus)
C
C**********************************************************************
C
C      Subroutine CVT_FNAME_DATA
C
C      Author:    Linda Wharton 1/93
C
C      Inputed DIR, GTIME and EXT are converted to ASCII values in a byte
C      array.  B_FNAME is created in C function make_c_fname.  FILE_NAME
C      is also created.  B_EXT is returned for use in other functions.
C
C      Modified:  Linda Wharton 4/94
C      Changed to remove references to BYTE data type, FILE_NAME is
C      created from DIR, GTIME and EXT_IN.  No ASCII conversions done.
C
C      Modified:  Linda Wharton 9/95
C      Removed downcase on whole filename...replaced by downcase on extension.
C
C      Modified:  Linda Wharton 2/96
C      Added second time variable, fhh, for making extended filenames
C      fhh is 4 characters, containing HHMM past rtime of validtime.
C      gtime variable changed to rtime for clarity.  Set up so only
C      when called with ext = LGA/lga will long filename be created.
C
C
C**********************************************************************

      IMPLICIT  NONE

      INTEGER*4 end_dir,end_ext,
     1          ERROR(2),
     1          ISTATUS

      INTEGER*2 fn_length

      CHARACTER*(*)      DIR                     !Directory to read data from
      CHARACTER*31      EXT             !File name ext (up to 31 chars)
      CHARACTER*31      EXT_IN          !File name ext (up to 31 chars) for downcase
      CHARACTER*9       RTIME
      CHARACTER*4       fhh
      CHARACTER*91      FILE_NAME

      ERROR(1)=1
      ERROR(2)=0

C **** Convert string data so it can be used by C programs
C ******  find end of dir and add to b_filname
C
      call s_len(dir,end_dir)
C
C ******  find end of ext
C
      call s_len(ext,end_ext)
C
C ****  make fortran file_name
C
      call downcase(ext,ext_in)

      if (ext_in .eq. 'lga') then
        file_name=dir(1:end_dir)//rtime//fhh//'.'//ext_in(1:end_ext)
        fn_length = end_dir + 14 + end_ext
      else
        file_name = dir(1:end_dir)//rtime//'.'//ext_in(1:end_ext)
        fn_length = end_dir + 10 + end_ext
      endif
C
C ****  Return normally.
C
        istatus=error(1)
999     return
        end

C########################################################################
      subroutine make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)

      implicit none

      integer*4       valtime, reftime, istatus
      integer*4       fcst_hr, fcst_min, fcst_min_sec, fcst_sec
      integer*4       error(3)
      character*4     fcst_hh_mm
      character*1       h1, h2, m1, m2

      error(1)=1
      error(2)=0

      fcst_sec = valtime - reftime
      if (fcst_sec .eq. 0) then
        fcst_hh_mm = '0000'
        goto 998
      endif

C ****  fcst_sec > 0 .... create fcst_hh_mm

      fcst_min_sec = mod(fcst_sec,3600)
      fcst_hr = (fcst_sec - fcst_min_sec) / 3600
      fcst_min = fcst_min_sec / 60

C     fcst_hr can be between 0 and 99
C     fcst_min can be between 0 and 59

      if ((fcst_hr .lt. 0) .or. (fcst_hr .gt. 99)) then
        write(6,*) ' Forecast hour in error: ',fcst_hr
        goto 997
      else
        if ((fcst_min .lt. 0) .or. (fcst_min .gt. 59)) then
          write(6,*) ' Forecast minute in error: ',fcst_min
          goto 997
        endif
      endif

      if (fcst_hr .lt. 10) then
        h1 = '0'
        h2 = char(fcst_hr + 48)
      else
        h1 = char(((fcst_hr - mod(fcst_hr,10)) / 10) + 48)
        h2 = char(mod(fcst_hr,10) + 48)
      endif

      if (fcst_min .lt. 10) then
        m1 = '0'
        m2 = char(fcst_min + 48)
      else
        m1 = char(((fcst_min - mod(fcst_min,10)) / 10) + 48)
        m2 = char(mod(fcst_min,10) + 48)
      endif

      fcst_hh_mm = h1//h2//m1//m2
      goto 998
C
C ****  Error Return.
C

997   istatus=error(2)
      goto 999

C
C ****  Return normally.
C

998   istatus=error(1)

999   return
      end
C########################################################################
