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
       subroutine read_laps(reftime,valtime,dir,ext,iimax,jjmax,
     1                      kkmax,kdim,var_req,lvl_req,lvl_coord_req,
     1                      units_req,comment_req,data,istatus)
C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            read_laps
C
C      The read_laps routine reads the following FORTRAN
C      routines from the readlapsdata.f file:
C            make_fcst_time
C            cvt_fname_v3
C
C      The read_laps routine calls the following C function
C      from the rwl_v3.c file:
C            read_cdf_v3
C
C**********************************************************************
C
C      Author:    Linda Wharton
C      Modified:  To accept netCDF ver. 2 data files  1/93 Linda Wharton
C                 To accept valtime & reftime  2/96 Linda Wharton
C                 To accept netCDF ver. 3 data files  9/97 Linda Wharton 
C
cdoc   Reads data requested by arrays VAR_REQ and LVL_REQ for the
cdoc   I4TIME, DIR and EXT specified.  Returns LVL_COORD-REQ,
cdoc   UNITS_REQ, COMMENT_REQ, DATA AND ISTATUS
cdoc   REFTIME is the time of the model run.
cdoc   VALTIME is the valid time of the data.
cdoc   For analysis, valtime = reftime.
C
C**********************************************************************
C
      implicit  none
C
      integer       reftime,             !INPUT I4time of model run
     1                valtime,             !INPUT I4time data is valid
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
      call make_fnam_lp(reftime,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'Error converting reftime to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)

      call s_len(ext, ext_len)

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930


      i_reftime = reftime - 315619200
      i_valtime = valtime - 315619200
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
     1    write (6,*) 'Error opening netCDF file...read aborted...'
     1               ,file_name
        istatus=error(2)
        goto 999
C
970     if (flag .ne. 1)
     1    write (6,*) 'Error retrieving data...read aborted.'
        istatus=error(2)
        goto 999
C
980     if (flag .ne. 1)
     1    write (6,*) 'Error in array dimensioning...read aborted.'
        istatus=error(2)
        goto 999
C
990     if (flag .ne. 1)
     1    write (6,*) 'Error in version, not a valid LAPS file...
     1read aborted.'
        istatus=error(2)
        goto 999
C
992     continue
        istatus=error(2)
        goto 999
C
        END

C########################################################################

