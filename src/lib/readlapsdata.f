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
       subroutine read_laps_data(i4time,dir,ext,iimax,jjmax,kkmax,kdim,
     1                    var_req,lvl_req,lvl_coord_req,
     1                    units_req,comment_req,data,istatus)
C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            read_laps_data
C            cvt_fname_v3
C            make_fcst_time
C
C      The read_laps_data subroutine calls the following C function
C      from the rwl_v3.c file:
C            read_cdf_v3
C
C**********************************************************************
C
C      Subroutine READ_LAPS_DATA
C
C      Author:    Linda Wharton
C      Modified:  To accept netCDF ver. 2 data files  1/93 Linda Wharton
C                 To accept netCDF ver. 3 data files...capablility to 
C                   read binary laps files separated out into 
C                   read_old_laps.f  9/97 Linda Wharton
C
C      Reads data requested by arrays VAR_REQ and LVL_REQ for the
C      I4TIME, DIR and EXT specified.  Returns LVL_COORD-REQ,
C      UNITS_REQ, COMMENT_REQ, DATA AND ISTATUS
C
C**********************************************************************
C
      implicit  none
C
      integer       i4time,              !INPUT I4time of data
     1                iimax,jjmax,kkmax,   !INPUT # cols, # rows, # fields
     1                kdim,                !INPUT K dimension of DATA array
     1                lvl_req(kdim),       !INPUT Requested levels
     1                istatus              !OUTPUT
      character*(*)   dir                  !INPUT Directory to read data from
      character*(*)   ext                  !INPUT File name ext 
      character*(*)   var_req(kdim)        !INPUT 3 letter ID of requested fields
      character*(*)   lvl_coord_req(kdim)  !OUTPUT Vertical coordinate of fields
      character*(*)   units_req(kdim)      !OUTPUT Units of requested fields
      character*(*)   comment_req(kdim)    !OUTPUT Comments for requested fields
      real        data(iimax,jjmax,kdim) !OUTPUT data
C
      integer fn_length,
     1          i_reftime,              !UNIX time of data
     1          i_valtime,              !UNIX time of data
     1          flag,                   !Print flag (1 = off)
     1          error(3),
     1          called_from,
     1          var_len,
     1          comm_len,
     1          ext_len,
     1          lvl_coord_len,
     1          units_len
  
C
      character*5       fcst_hh_mm
      character*9       gtime
      character*150     file_name
C
      common            /prt/flag
C
C-------------------------------------------------------------------------------
C
      error(1)=1
      error(2)=0
      error(3)=-2
C
C ****  Create file name.
C
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'Error converting i4time to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call s_len(ext, ext_len)

C Hard wired as a place holder - will be used in filename only if read_laps_data
C   is called on lga, lgb, fua, fsf, ram, rsf
C To fix this, call read_laps instead
      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930
C
C **** get actual reftime from gtime...
C
      i_reftime = i4time - 315619200
      i_valtime = i_reftime

      called_from = 0   !called from FORTRAN

      var_len = len(var_req(1))
      comm_len = len(comment_req(1))
      lvl_coord_len = len(lvl_coord_req(1))
      units_len = len(units_req(1))

      call read_cdf_v3 (file_name, ext, var_req, comment_req, 
     1                  lvl_coord_req, units_req, var_len, comm_len, 
     1                  fn_length, ext_len, lvl_coord_len, units_len, 
     1                  i_reftime, i_valtime, iimax, jjmax, kkmax,  
     1                  kdim, lvl_req, data, called_from, istatus)


      if (istatus .ge. 0) then   !return from read with no errors
                                 !convert byte data to characters

         if (istatus .gt. 0) then
            istatus=error(3)
         else
            istatus=error(1)
         endif
         goto 999
      endif

      if (istatus .eq. -5)  goto 992  !file not there
      if (istatus .eq. -4)  goto 990  !error in version 
      if (istatus .eq. -3)  goto 980  !error in dimensioning arrays
      if (istatus .eq. -2)  goto 970  !error retrieving data
      if (istatus .eq. -1)  goto 950  !error opening file as netCDF
C
C ****  Return normally.
C
        istatus=error(1)
999     return
C
C ****  Error trapping.
C
930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...read aborted.'
        istatus=error(2)
        goto 999
C
950     if (flag .ne. 1)
     1    write (6,*) 'Error opening netCDF file...read aborted.'
        istatus=error(2)
        goto 999
C
970     if (flag .ne. 1)
     1    write (6,*) 'Error retrieving data...read aborted.'
        istatus=error(2)
        goto 999
C
980     if (flag .ne. 1) then
          write (6,*) 'Error in array dimensioning...read aborted.'
          write (6,*) 'iimax/jjmax/kkmax/kdim = ',iimax,jjmax,kkmax,kdim
        endif
        istatus=error(2)
        goto 999
C
990     if (flag .ne. 1)
     1    write (6,*) 'Error in version, not a valid LAPS file... '
     1                ,'read aborted.'
        istatus=error(2)
        goto 999
C
992     continue
        istatus=error(2)
        goto 999
C
        END

C########################################################################
      subroutine cvt_fname_v3(dir,gtime,fhh,ext,ext_len,file_name,
     1                        fn_length,istatus)
C**********************************************************************
C
C      Subroutine CVT_FNAME_V3
C
C      Author:    Linda Wharton 1/93
C
C      Inputed DIR, GTIME and EXT are converted to ASCII values in a byte
C      array.  B_FNAME is created in C function make_c_fname.  FILE_NAME
C      is also created.  B_EXT is returned for use in other functions.
C
C      Modified:  Linda Wharton 4/94
C      Changed to remove references to BYTE data type, FILE_NAME is
C      created from DIR, GTIME and EXT.  No ASCII conversions done.
C
C**********************************************************************

      implicit  none

      character*(*)     dir          !Directory to read data from
      character*(*)     gtime
      character*(*)     ext          !File name ext 
      character*(*)     file_name
      character*5       fhh

      integer         fn_length,
     1                  ext_len,fhh_len,
     1                  istatus

      integer         end_dir, end_ext, error(2)

C#ifdefined NODYNAMIC
C      character*31  ext_dn
C#else
      character*(10) ext_dn
C#endif

      call downcase(ext,ext_dn)

      error(1)=1
      error(2)=0

C **** Convert string data so it can be used by C programs
C ******  find end of dir 
C
      call s_len(dir,end_dir)
C
C ******  find end of ext_dn
C
      call s_len(ext_dn,end_ext)
C
C ******  find end of file_name
C
      fn_length = len(file_name)
C
C ******  find end of fhh
C
      call s_len(fhh,fhh_len)
C
C ****  make fortran file_name
C

      if (end_dir+end_ext+fhh_len+10 .gt. fn_length) then
        write (6,*) 'Length of dir+file-name exceeds file_name length.'
        istatus = error(2)
        goto 999
      else
        if (ext_dn(1:2) .eq. 'lg' .or.
     +     ext_dn(1:3).eq.'fua' .or.
     +     ext_dn(1:4).eq.'cont' .or.
     +     ext_dn(1:3).eq.'fsf' .or.
     +     ext_dn(1:3).eq.'ram' .or.
     +     ext_dn(1:3).eq.'rsf') then
           file_name=dir(1:end_dir)//gtime//fhh(1:fhh_len)//'.'//
     +ext_dn(1:end_ext)
        else
          file_name = dir(1:end_dir)//gtime//'.'//ext_dn(1:end_ext)
        endif

        call s_len(file_name,fn_length)
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

      integer       valtime, reftime, istatus
      integer       fcst_hr, fcst_min, fcst_min_sec, fcst_sec
      integer       error(3), extended, interim_fh
      character*5   fcst_hh_mm
      character*1   h1, h2, h3, m1, m2

      error(1)=1
      error(2)=0

      fcst_sec = valtime - reftime

C see if fcst_sec > 356400 seconds (99 hours)
      if (fcst_sec .gt. 356400) then
        extended = 1  ! fcst_hh_mm format is hhhmm if hh > 99
      else
        extended = 0  ! fcst_hh_mm format is hhmm if hh <= 99
      endif

      if (fcst_sec .eq. 0) then
        fcst_hh_mm = '0000'
        goto 998
      endif

C ****  fcst_sec > 0 .... create fcst_hh_mm

      fcst_min_sec = mod(fcst_sec,3600)
      fcst_hr = (fcst_sec - fcst_min_sec) / 3600
      fcst_min = fcst_min_sec / 60

C     fcst_hr can be between 0 and 999
C     fcst_min can be between 0 and 59

      if ((fcst_hr .lt. 0) .or. (fcst_hr .gt. 999)) then
       write(6,*) ' Forecast hour cannot exceed 999: ',fcst_hr
        goto 997
      else
        if ((fcst_min .lt. 0) .or. (fcst_min .gt. 59)) then
          write(6,*) ' Forecast minute cannot exceed 59: ',fcst_min
          goto 997
        endif
      endif

      if (extended .eq. 1) then  !fcst_hr > 99, so 3 digit hr is valid
        interim_fh = fcst_hr - mod(fcst_hr,100) 
        h1 = char((interim_fh/100) + 48)
        interim_fh = fcst_hr - interim_fh 
        h2 = char(((interim_fh - mod(interim_fh,10)) / 10) + 48)
        h3 = char(mod(interim_fh,10) + 48)
      else ! fcst_hr <= 99, so 2 digit hr is valid
        if (fcst_hr .lt. 10) then
          h1 = '0'
          h2 = char(fcst_hr + 48)
        else
          h1 = char(((fcst_hr - mod(fcst_hr,10)) / 10) + 48)
          h2 = char(mod(fcst_hr,10) + 48)
        endif
      endif

      if (fcst_min .lt. 10) then
        m1 = '0'
        m2 = char(fcst_min + 48)
      else
        m1 = char(((fcst_min - mod(fcst_min,10)) / 10) + 48)
        m2 = char(mod(fcst_min,10) + 48)
      endif

      if (extended .eq. 1) then
        fcst_hh_mm = h1//h2//h3//m1//m2
      else
        fcst_hh_mm = h1//h2//m1//m2
      endif
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
