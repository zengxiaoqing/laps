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
      SUBROUTINE READ_LAPS_HEADER(I4TIME,DIR,EXT,IMAX,JMAX,KMAX,
     1                         LAPS_DOM_FILE,ASCTIME,VERSION,
     1                         MODEL,ORIGIN,VAR,LVL,NUM_VARIABLES,
     1                         VAR_AVAIL,LAPS_VAR_AVAIL,NUM_LEVELS,
     1                         LVL_AVAIL,LVL_COORD,UNITS,
     1                         COMMENT,L_PACKED_DATA,ISTATUS)
C
C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            readlapsheader
C
C      The readlapsheader subroutine reads the following FORTRAN
C       subroutines from the readlapsdata.for file:
C            cvt_fname_data
C
C      The readlapsheader subroutine reads the following C subroutines
C      from the readwritelaps.c file:
C            make_c_fname
C            read_cdf_header
C            cdf_retrieve_header
C            itoa
C
C**********************************************************************
C
C      Subroutine READ_LAPS_HEADER
C
C      Author:    Steve Albers
C      Modified:  To accept netCDF data files          1/93 Linda Wharton
C                 To accept extended format filenames  2/96  Linda Wharton
C                 Will not make correct filetimes with
C                 extension lga
C
C      Creates filename using cvt_fname_data.
C      Reads header variables IMAX, JMAX, KMAX, LAPS_DOM_FILE, ASCTIME,
C      VERSION, MODEL, ORIGIN, and NUM_VARIABLES.
C
C**********************************************************************
        IMPLICIT        NONE

        INTEGER         FN_LENGTH
C
        INTEGER       I4TIME,         ! Input I4time of data
     1          FLAG,           ! Print flag (1 = off)
     1          no_laps_diag,   !if = 0, print diagnostic output
     1          REC,
     1          RCDL,
     1          IMAX,JMAX,KMAX, ! Output dimensions/#fields from header
     1          LVL(200),      ! Output levels from header
     1          I,K,
     1          ERROR(2),
     1                num_variables,
     1                NUM_LEVELS,
     1          ISTATUS
C
        REAL          MSG_FLAG
C
        CHARACTER*150   DIR             ! Input Directory to read data from
        CHARACTER*31    EXT             ! Input File name ext (up to 31 chars)
        CHARACTER*31    EXT_I           !INPUT input file name ext
        CHARACTER*3     VAR(200)       ! Output Variables
        CHARACTER*4     LVL_COORD(200) ! Output comments
        CHARACTER*10    UNITS(200)     ! Output units
        CHARACTER*125   COMMENT(200)   ! Output comments
        CHARACTER*9     GTIME
        CHARACTER*5     fcst_hh_mm
        CHARACTER*91    FILE_NAME
        CHARACTER*4     MARK
        CHARACTER*4     CIMAX,CJMAX,CKMAX
        CHARACTER*24    ASCTIME
        character*18    asct
        CHARACTER*4     VERSION
        character*5     vern
        CHARACTER*131   MODEL           !Meteorological model in file
        character*132   modl
        CHARACTER*131   ORIGIN          !Location where file was created
        character*132   orign
        CHARACTER*11    LAPS_DOM_FILE   !Name of domain file e.g. NEST7GRID
        character*12    ldf
        CHARACTER*1     HMARK
        CHARACTER*4     CLVL
        CHARACTER*4     CSTART_REC
        CHARACTER*3     LAPS_VAR_AVAIL(200)
        character*4     lvar_a(200)
        CHARACTER*19    VAR_AVAIL(200)
        character*20    var_a(200)
C
        LOGICAL         l_packed_data
C
        INTEGER       LVL_AVAIL(200)
C
        DATA            MSG_FLAG/9.E30/
C
        COMMON          /PRT/FLAG
        COMMON          /laps_diag/no_laps_diag

C
C-------------------------------------------------------------------------------
C
        ERROR(1)=1
        ERROR(2)=0
C
C ****  Create file name.
C
        CALL MAKE_FNAM_LP(I4TIME,GTIME,ISTATUS)
        IF (ISTATUS .ne. 1) THEN
                ISTATUS=ERROR(2)
                RETURN
        ENDIF

C fcst_hh_mm: Hard wired as a place holder - will be used in filename 
C   only if read_laps_header is called on lga, lgb, fua, fsf, ram, rsf
        fcst_hh_mm = '0000'

        CALL UPCASE(EXT,EXT_I)

        call cvt_fname_data(dir,gtime,fcst_hh_mm,ext_i,file_name,
     1                      fn_length,istatus)

        call read_cdf_header(file_name, fn_length,imax, jmax, kmax,
     1                 num_variables,ldf,asct,vern,modl,orign,
     1                 var_a,lvar_a,num_levels, lvl_avail,
     1                 no_laps_diag,Istatus)

        if (istatus .eq. 1) then
          laps_dom_file = ldf
          asctime = asct
          version = vern
          model = modl
          origin = orign
          do i = 1, num_variables
            laps_var_avail(i) = lvar_a(i)
            var_avail(i) = var_a(i)
          enddo

         ISTATUS=ERROR(1)
         call setup_var_lvl(ext_i,num_levels,lvl_avail,
     1                      num_variables,LAPS_var_avail,
     1                      var,lvl,kmax,istatus)
         if (istatus .ne. 1) then
            write (6,*) 'Error in setup_var_lvl'
            ISTATUS=ERROR(2)
         endif
         l_packed_data = .false.
         goto 999
        else
          if (istatus .eq. -2) goto 905

!         if you get to here, istatus = -1 so file not netCDF

c         IF (ext(1:3) .eq. 'LC3' .or. ext(1:3) .eq. 'SC3' .or.
c    1    ext(1:3) .eq. 'LPS' .or. ext(1:3) .eq. 'SPS' .or.
c    1    ext(1:3) .eq. 'LRP' .or. ext(1:3) .eq. 'SRP' .or.
c    1    ext(1:3) .eq. 'LTY' .or. ext(1:3) .eq. 'STY' .or.
c    1    ext(1:3) .eq. 'LH3' .or. ext(1:3) .eq. 'SH3' .or.
c    1    ext(1:3) .eq. 'LMD' .or. ext(1:3) .eq. 'SMD'            )then
            l_packed_data = .false.
c         else
c           l_packed_data = .false.
c         endif
          if(l_packed_data)then
             goto 905
          endif
C
C ****  Open file and read first header record.
C
        OPEN(1,FILE=FILE_NAME,STATUS='OLD',ERR=940)
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
        IF(KMAX .GT. 200)GOTO 960
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
        ENDDO
906     FORMAT(I3)

910     FORMAT(1X,A3,I4,' field not found.')

        endif  !if status .eq. 1
C
C ****  Return normally.
C
        l_packed_data = .false.
        ISTATUS=ERROR(1)
998     CLOSE(1,ERR=999)
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
     1    write (6,*) 'KMAX greater than 200...arrays to small'
        ISTATUS=ERROR(2)
        GOTO 998
C
        END

C########################################################################
      subroutine setup_var_lvl(ext_in,num_levels,lvl_avail,
     1                         num_variables,LAPS_var_avail,
     1                         var,lvl,kdim,istatus)

C**********************************************************************
C
C      Subroutine SETUP_VAR_LVL
C
C      Author:    Linda Wharton
C
C      Takes data read in by READLAPSHEADER and fills VAR and LVL
C      so READLAPSDATA may be called.  This subroutines is used only
C      if the file is in netCDF format.
C
C**********************************************************************
C
      IMPLICIT  NONE

      INTEGER LVL(*),
     1          i, j,
     1          KDIM,
     1          NUM_LEVELS,
     1          num_variables,
     1          ERROR(3),
     1          pos,
     1          ISTATUS                     !OUTPUT

      INTEGER  LVL_AVAIL(*)

      CHARACTER*3       LAPS_VAR_AVAIL(*)
      CHARACTER*3       VAR(*)
      CHARACTER*31      EXT_IN              !INPUT input file name ext

C
C-------------------------------------------------------------------------------
C
      ERROR(1)=1
      ERROR(2)=0
      ERROR(3)=-2

      istatus = error(1)

      if (ext_in .eq. 'LMR') then
         var(1) = 'R00'
         var(2) = 'R06'
         var(3) = 'R12'
         lvl(1) = 0
         lvl(2) = 0
         lvl(3) = 0
         pos = 4
      else
         if (ext_in .eq. 'LF1') then
            var(1) = 'H00'
            var(2) = 'H06'
            var(3) = 'H12'
            lvl(1) = 0
            lvl(2) = 0
            lvl(3) = 0
            pos = 4
          else
           if (ext_in .eq. 'LHE') then
              var(1) = 'LHE'
              lvl(1) = 0
              pos = 2
           else
              pos = 1

              do i = 1, num_variables
                 if (LAPS_var_avail(i)(1:1) .ne. char(0)) then
                    do j = 1, num_levels
                       var(pos) = LAPS_var_avail(i)
                       lvl(pos) = lvl_avail(j)
                       pos = pos + 1
                    enddo
                 endif
              enddo
            endif
         endif
      endif

      if ((pos - 1) .ne. kdim) istatus = error(2)

      return
      end

C########################################################################
