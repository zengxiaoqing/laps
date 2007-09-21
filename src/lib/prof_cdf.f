cdis
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
C       PROF_CDF_CLOSE.FOR              Michael Barth           11-Aug-1993
C
C       This subroutine is used to close a netCDF file of profiler data.
C
C INPUT:
C
C       CDFID           INTEGER:  The netCDF id of the file to be closed.  It
C                       must previously have been opened by PROF_CDF_OPEN.
C
C       STATUS          INTEGER:  Returned status:  0 is good, -1 and +n are
C                       errors returned by netCDF.  The definitions of these can
C                       be found in NETCDF.INC.  In addition, these error
C                       returns are from PROF_CDF_CLOSE:
C
C                       -5      NetCDF file not open.
C
C
C Modifications:
C       28-Nov-1995/LAB Put prof_cdf_common in it's own include file.
C       07-Aug-2000/LSW modified to netCDF version 3.4

C
        subroutine PROF_CDF_CLOSE(cdfid,status)
C
        implicit none
C
        integer cdfid,status
C
C       Internal PROF_CDF information:
C
C       OPEN_LIST       Logical flags indicating which of the four file slots
C                       are in use.
C
C       CDFID_LIST      NetCDF id's for the close files.
C
C       STANAM_LIST     Station names for the records in the close files.
C
C       WMOID_LIST      WMO id's for the records in the close files.
C
C       NREC_LIST       Number of records (stations) in each close file.
C
C       ERROR_FLAG      Error handling flag.
C
        integer max_profilers
        parameter (max_profilers = 200)
C
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
C
        integer i
C
        do i = 1, 4
C
C       Find the user's CDFID in our internal list.
C
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
C
                status = NF_CLOSE(cdfid_list(i))
C
                open_list(i) = .false.                  ! Indicate internally.
                go to 900
           endif
        enddo
C
        status = -5                                     ! File not found.
C
C       Use the non-default error handling for PROF_CDF errors.
C
        if(error_flag.ne.0)then
                write(*,*)'PROF_CDF_CLOSE:  NetCDF file not open.'
                if(error_flag.eq.2)stop
        endif
C
900     return
        end

C   PROF_CDF_OPEN.FOR           Michael Barth           11-Aug-1993
C
C       This subroutine is used to open a netCDF file holding profiler data.
C       It will also read in the station names and WMO ID's for later use when 
C       read routines attempt to locate data.
C
C       A user of the PROF_CDF routines is allowed to have up to 4 open
C       netCDF files at a time.
C
C INPUT:
C
C       FILE            Character string giving the name (including any path
C                       information) of the netCDF file to be opened.
C
C OUTPUT:
C
C       CDFID           INTEGER:  The netCDF id for the open file.
C
C       STATUS          INTEGER:  Returned status:  0 is good, -1 and +n are
C                       errors returned by netCDF.  The definitions of these can
C                       be found in NETCDF.INC.  In addition, these error
C                       returns are from PROF_CDF_OPEN:
C
C                       -2      Attempt to open more than 4 files.
C
C MODIFICATIONS:
C
C       28-Nov-1995/LAB Put prof_cdf_common in it's own include file.
C
C       07-Mar-1996/MFB Added the use of RECDIM_LIST to determine if a
C                       one-dimensional variable's dimension is the record
C                       dimension or not -- the hyperslab COUNT
C                       depends on this distinction. ***** Note that the
C                       RECDIM_LIST calculation is currently set to
C                       compensate for a bug in the old VMS netCDF library
C                       in use.  See below where RECDIM_LIST is calculated
C                       when porting this routine to Unix (or in changing
C                       netCDF library versions).
C
C       21-Mar-1996/MFB Finally went out and got the latest version of
C                       netCDF (2.4.1) and installed on my cluster.  This
C                       eliminated the need to compensate for the netCDF
C                       bug described in the 3/7/96 modification note.
C                       Also removed declarations that duplicate those in new
C                       netcdf.inc.
C
C       01-May-1996/MFB Went back to the old version of the netCDF library
C                       due to performance problems with the new one.  Thus
C                       the bug in NCINQ mentioned in the 3/7/96 mod note
C                       again has to be accommodated.  See the code below
C                       for more explanation.
C
C       26-Feb-1997/MFB Fixed a bug that was exposed if an existing file
C                       is opened, but it has no records in it.  This should
C                       be legal, and now that I've tried it and found the
C                       bug, it is legal...
C
C       28-Feb-1997/LAB Removed a forgotten write statement
C
C       17-Mar-1997/MFB Removed the RECDIM fix for VMS's old netCDF library.
C                       This should fix problems currently encountered with
C                       single dimensional variables.
C
C       16-Apr-1997/MFB Removed duplicate definitions of netCDF functions
C                       that are contained in netcdf.inc.
C
C       27-May-1997/MFB Renamed "ncopts" common and variable so it doesn't
C                       conflict with symbols used in C netCDF interface.
C
C       07-Aug-2000/LSW modified to netCDF version 3.4
C
        subroutine PROF_CDF_OPEN(file,cdfid,status)
C
        implicit none
C
        character*(*) file
        integer cdfid,status
C
C       Internal PROF_CDF information:
C
C       OPEN_LIST       Logical flags indicating which of the four file slots
C                       are in use.
C
C       CDFID_LIST      NetCDF id's for the open files.
C
C       STANAME_LIST    Station names for the records in the open files.
C
C       WMOID_LIST      WMO id's for the records in the open files.
C
C       NREC_LIST       Number of records (stations) in each open file.
C
C       RECDIM_LIST     Record dimension ID's for each file.
C
C       ERROR_FLAG      Error handling flag.
C
        integer max_profilers
        parameter (max_profilers = 200)
C
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
        character*128 dimname                   ! Must match NETCDF.INC's
C                                               ! MAXNCNAM parameter.
        integer ncopts_val
        common/ncopts_cmn/ncopts_val            ! NetCDF error handling flag.
C
        integer i,j,varid,start(2),count(2)
        integer ndims,nvars,ngatts,recdim

        logical l_exist

        write(6,*)' Subroutine PROF_CDF_OPEN'

        inquire(file=file,exist=l_exist)
        if(.not. l_exist)then
            write(6,*)' Warning in PROF_CDF_OPEN, file does not exist: '
     1               ,file 
            status=-7
            return
        endif ! l_exist
C
C       Unless non-default error handling has been set by a call to
C       PROF_CDF_SET_ERROR set the default error handling.
C
C
        do i = 1, 4
C
C       Find an open slot.
C
           if(.not.open_list(i))then
C
C       Open the file.
C
              status = NF_OPEN(file,NF_NOWRITE,cdfid)
              if(status.eq.NF_NOERR)then
C
C       Clear out the data structures for this slot.
C
                 do j = 1, max_profilers
                    staname_list(j,i) = '      '
                    wmoid_list(j,i) = 0
                 enddo
                 nrec_list(i) = 0
C
C       Get the number of records (stations) in the file.
C
                 status = NF_INQ_DIMID(cdfid, 'recNum',varid)
                 if (status.eq.NF_NOERR)then
                   status = NF_INQ_DIMLEN(cdfid,varid,nrec_list(i))

                 
                   if(nrec_list(i) .gt. max_profilers)then
                     write(6,*)' Too many profilers in PROF_CDF_OPEN '       
     1                        ,nrec_list(i),max_profilers
                     status=-7
                     return
                   endif
                 endif

                 if(status.eq.NF_NOERR.and.nrec_list(i).gt.0)then
C
C       Get the station names.  (Have to do this in a loop as I can't get
C       NCVGTC to work on one call for the whole enchilada.)
C
                    status = NF_INQ_VARID(cdfid,'staName',varid)

                    start(1) = 1
                    count(1) = 6 ! 6 characters
                    count(2) = 1 ! 1 record at a time
                    do j = 1, nrec_list(i)
                       start(2) = j
                       status = NF_GET_VARA_TEXT(cdfid,varid,start,
     $                                           count,
     $                                           staname_list(j,i))
                    enddo
                    if(status.eq.0)then
C
C       Get the WMO id's.
C
                       start(1) = 1 ! Start at the beginning
                       start(2) = 1
                       count(1) = nrec_list(i) ! and get all records.
                       status = NF_INQ_VARID(cdfid,'wmoStaNum',varid)
                       status = NF_GET_VARA_INT(cdfid,varid,start,count,
     $                                          wmoid_list(1,i))
                       if(status.ne.NF_NOERR)go to 900
                    endif
                 endif
C
C       Get the dimension id of the record dimension into the list where it
C       can be used by PROF_CDF_READ.  Note the +1 --> in the old netCDF
C       bug-riddled library NCINQ returns 0-(NDIMS-1) but the FORTRAN interface
C       requires 1-NDIM.  Note that in the current version of the netCDF
C       library (4/96) this bug has been fixed, but other problems exist.  Thus
C       we're sticking with the old version until porting this to Unix.
C
C       3/17/97: Have ported to Unix, and have removed the +1.  If this rears
C       it's head again the OLD code was: recdim_list(i) = recdim + 1
C
                status = NF_INQ(cdfid,ndims,nvars,ngatts,recdim)
                if(status.eq.NF_NOERR)recdim_list(i) = recdim
C
              endif
              go to 900                 ! Done processing with the open slot.
           endif
        enddo
C
        status = -2                     ! No open slots.
C
900     if(status.eq.0)then
C
           open_list(i) = .true.
           cdfid_list(i) = cdfid
C
        else
C
           cdfid = 0
C
           if(error_flag.gt.0.and.status.eq.-2)then
C
              write(*,*)'PROF_CDF_OPEN:  Attempt to open more than 4'
     $                  //' files.'
C
           endif
C
 
           if(error_flag.eq.2)stop
 
C
        endif
C
        return
        end
C       PROF_CDF_READ.FOR               Michael Barth           11-Aug-1993
C
C       This subroutine is used to read all of one variable for one station
C       from a netCDF file holding profiler data.
C
C       The user must have previously opened the file via PROF_CDF_OPEN.
C
C INPUT:
C
C       CDFID           INTEGER:  The netCDF id for the file.
C
C       STANAME         CHARACTER*6:  The name of the station (e.g., 'PLTC2 ' --
C                       note the trailing blank).  This is the primary way the
C                       user can request a station.  If this is all blanks,
C                       the second method, WMOID, will be used.
C
C                       ***** Note that the station name must be all capital
C                       letters (and a number).
C
C       WMOID           INTEGER:  The WMO block and station number (e.g.,
C                       74533).  This is the alternative method of selecting
C                       the station.
C
C                       ***** Note that if STANAME is not all blank, WMOID
C                       will be ignored.
C
C       VARNAME         CHARACTER*(*):  A character string giving the name of
C                       the desired variable (use the names in the netCDL
C                       file describing the dataset in use).  This is
C                       case-sensitive so the exact string must be used.
C
C       DTYPE           INTEGER:  0 - The variable is numeric data.
C                                   1 - The variable is character data.
C
C OUTPUT:
C
C       ARRAY           Array of data returned to the caller.  The data type
C                       and dimensions are not used in this routine, or
C                       validated.  The caller must use the appropriate type
C                       (see the netCDL file) and know the expected amount
C                       of data.
C
C       STATUS          INTEGER:  Returned status:  0 is good, -1 and +n
C                       are errors returned by netCDF.  The definitions of
C                       these can be found in NETCDF.INC.  In addition, these
C                       error returns are from PROF_CDF_READ:
C
C                       -3      Station not found.
C                       -4      Invalid DTYPE.
C                       -5      NetCDF file not open.
C                       -9      Too many dimensions in variable (16 maximum).
C
C MODIFICATIONS:
C
C       15-Sep-1994/MFB Fixed a really stupid oversight:  the routine didn't
C                       work with 2-dimensional arrays for a single record
C                       (e.g., moments - rec_num, level, beam).  Part of this
C                       fix was to remove the NARRAY argument into this
C                       routine.  That was the size of the expected array.
C                       As the fix included the need to determine that on a
C                       dimension-by-dimension basis, NARRAY became redundant.
C
C       11-Sep-1995/MFB Change the dimensions of START & COUNT to use maximum
C                       number of allowable netCDF dimensions.  Also changed
C                       how character data are read to fix a bug revealed by
C                       the new BLP CDL.
C
C       28-Nov-1995/LAB Put prof_cdf_common in it's own include file.
C
C       29-Feb-1996/MFB Changed declaration of ARRAY from BYTE to INTEGER to
C                       make the code more portable.
C
C       05-Mar-1996/MFB Limited the check of STANAME to the first five
C                       characters.  This was done to make the code compatible
C                       with all known profiler-based data formats in FSL.
C                       (Some have a blank in the 6th character, some have a
C                       zero.)
C
C       07-Mar-1996/MFB Added the use of RECDIM_LIST to determine if a
C                       one-dimensional variable's dimension is the record
C                       dimension or not -- the hyperslab COUNT
C                       depends on this distinction.
C
C       21-Mar-1996/MFB Removed declarations that duplicate those in new
C                       netcdf.inc.
C
C       07-Aug-2000/LSW modified to netCDF version 3.4
C
        subroutine PROF_CDF_READ(cdfid,staname,wmoid,varname,dtype,
     $                           array,status)
C
        implicit none
C
        integer cdfid,wmoid,dtype,status
        character*6 staname
        character*(*) varname
        integer array(*)
C
C       Internal PROF_CDF information:
C
C       OPEN_LIST       Logical flags indicating which of the four file slots
C                       are in use.
C
C       CDFID_LIST      NetCDF id's for the read files.
C
C       STANAME_LIST    Station names for the records in the read files.
C
C       WMOID_LIST      WMO id's for the records in the read files.
C
C       NREC_LIST       Number of records (stations) in each read file.
C
C       RECDIM_LIST     Record dimension ID's for each file.
C
C       ERROR_FLAG      Error handling flag.
C
        integer max_profilers
        parameter (max_profilers = 200)
 
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
C
        character*21 error_text(3)
        data error_text/'Station not found.   ',
     $                  'Invalid DTYPE.       ',
     $                  'NetCDF file not open.'/
        character*(maxncnam) name
        integer vartyp,nvdims,vdims(maxvdims),nvatts,dimsiz
        integer i,j,k,varid,start(maxvdims),count(maxvdims)
        integer count_total
        character*5 c5_1,c5_2
C
        do i = 1, 4
C
C       Find the user's CDFID in our internal list.
C
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
C
C       Locate the station.
C
              do j = 1, nrec_list(i)
C
                 c5_1 = staname_list(j,i)(1:5)
                 c5_2 = staname(1:5)
C
                 if((staname.ne.'      '.and.c5_1.eq.c5_2).or.
     $              (staname.eq.'      '.and.
     $               wmoid_list(j,i).eq.wmoid))then
C
C       Get the variable id.
C
                    status = NF_INQ_VARID(cdfid,varname,varid)
                    if(status.ne.NF_NOERR)go to 900
C
C       Find out how many dimensions are in the variable.
C
                    status = NF_INQ_VAR(cdfid,varid,name,vartyp,
     $                                  nvdims,vdims,nvatts)
                    if(status.ne.NF_NOERR)go to 900
C
C       Use this information to, one dimension at a time, fill in the
C       dimensions of the hyperslab to be obtained.
C
                    count_total = 1     ! Must calculate product of COUNT
C                                       ! vector for character data type.
                    do k = 1, nvdims
                       status = NF_INQ_DIM(cdfid,vdims(k),name,dimsiz)
                       if(status.ne.NF_NOERR)go to 900
                       if(k.lt.nvdims.or.
     $                    (nvdims.eq.1.and.vdims(k).ne.
     $                     recdim_list(i)))then
                                start(k) = 1
                                count(k) = dimsiz
                                count_total = count_total * dimsiz
                       else
                                start(k) = j    ! REC_NUM dimension
                                count(k) = 1
                       endif
                    enddo
C
C       Call the appropriate routine to get the variable.
C
                    if(dtype.eq.0)then  !reading an integer
C
                      status = NF_GET_VARA_INT(cdfid,varid,start,count,
     $                                         array)
C
                    else if(dtype.eq.1)then   !reading text
C
                      status = NF_GET_VARA_TEXT(cdfid,varid,start,count,
     $                                          array)
C
                    else if(dtype.eq.2)then   !reading real/float
C
                      status = NF_GET_VARA_REAL(cdfid,varid,start,count,
     $                                          array)
                    else
C
                             status = -4                ! Invalid DTYPE.
C
                    endif
                    if (status.ne.NF_NOERR) print *, NF_STRERROR(status)
C
                    go to 900                   ! End of found station.
C
                 endif
              enddo
C
              status = -3                       ! Station not found.
              go to 900
C
           endif                                ! End of found open file.
C
        enddo
C
        status = -5                             ! File not found.
C
900     if(status.lt.-1.and.error_flag.ne.0)then
C
C       Use the non-default error handling for PROF_CDF errors.
C
           
                write(*,*)'PROF_CDF_READ:  '//
     $                    error_text(abs(status)-2)
                if(error_flag.eq.2)stop
        endif
C
        return
        end
C       PROF_CDF_SET_ERROR.FOR          Michael Barth           12-Aug-1993
C
C       This subroutine is used to set error handling for PROF_CDF routines.
C       The caller doesn't have to use this feature, if no call is made the
C       default error handling will be performed.
C
C INPUT:
C
C       ERROR_CODE      0       Return status codes -- THIS IS THE DEFAULT.
C                       1       Return status codes and write an error message
C                               to standard output (SYS$OUTPUT on VMS).
C                       2       Write an error message to standard output and
C                               exit the program.
C
C OUTPUT:
C
C       STATUS          INTEGER:  Returned status:  0 is good, -1 and +n are
C                       errors returned by netCDF.  The definitions of these can
C                       be found in NETCDF.INC.  In addition, these error
C                       returns are from PROF_CDF_SET_ERROR:
C
C                       -6      Invalid ERROR_CODE.
C
C MODIFICATIONS:
C
C       28-Nov-1995/LAB Put prof_cdf_common in it's own include file.
C
C       27-May-1997/MFB Renamed "ncopts" common and variable so it doesn't
C                       conflict with symbols used in C netCDF interface.
C
C       07-Aug-2000/LSW modified to netCDF version 3.4 NCPOPT is no longer
C                       supported, and as such is no longer called.  
C
        subroutine PROF_CDF_SET_ERROR(error_code,status)
C
        implicit none
C
        integer error_code,status
C
C       Internal PROF_CDF information:
C
C       OPEN_LIST       Logical flags indicating which of the four file slots
C                       are in use.
C
C       CDFID_LIST      NetCDF id's for the open files.
C
C       STANAM_LIST     Station names for the records in the open files.
C
C       WMOID_LIST      WMO id's for the records in the open files.
C
C       NREC_LIST       Number of records (stations) in each open file.
C
C       ERROR_FLAG      Error handling flag.
C
        integer max_profilers
        parameter (max_profilers = 200)
C
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
C
        integer ncopts_val
        common/ncopts_cmn/ncopts_val            ! NetCDF error handling flag.
C
C       Validate the ERROR_CODE and then save it in the PROF_CDF common.  Then
C       set the netCDF version (NCOPTS_VAL) accordingly.
C
        status = 0                              ! Assume success.
C
        if(error_code.eq.0)then
C
                error_flag = 0
 
C
        else if(error_code.eq.1)then
C
                error_flag = 1
C
        else if(error_code.eq.2)then
C
                error_flag = 2
        else
                status = -6                     ! Invalid ERROR_CODE.
        endif
C
        return
        end
C       PROF_CDF_GET_STATIONS.FOR       Michael Barth           19-Oct-1995
C
C       This subroutine is used to get a list of stations from a netCDF file
C       holding profiler data.
C
C       The user must have previously opened the file via PROF_CDF_OPEN.
C
C INPUT:
C
C       CDFID           INTEGER:  The netCDF id for the file.
C
C OUTPUT:
C
C       N_STATIONS      Number of stations currently in the file.
C
C       STANAME         CHARACTER*6 array:  The names of the stations (e.g.,
C                       'PLTC2 ').
C
C       WMOID           INTEGER array :  The WMO block and station numbers
C                       (e.g., 74533).
C
C       STATUS          INTEGER:  Returned status:  0 is good.  No netCDF
C                       calls are performed, so no errors are returned by
C                       netCDF.  In addition, these error returns are from
C                       PROF_CDF_GET_STATIONS:
C
C                       -5      NetCDF file not open.
C
C MODFICATIONS:
C       28-Nov-1995/LAB Put prof_cdf_common in it's own include file.
C
C       07-Aug-2000/LSW modified to netCDF version 3.4
C
        subroutine PROF_CDF_GET_STATIONS(cdfid,n_stations,staname,
     $                                   wmoid,status)
C
        implicit none
C
        integer cdfid,n_stations,wmoid(*),status
        character*6 staname(*)
C
C       Internal PROF_CDF information:
C
C       OPEN_LIST       Logical flags indicating which of the four file slots
C                       are in use.
C
C       CDFID_LIST      NetCDF id's for the read files.
C
C       STANAM_LIST     Station names for the records in the read files.
C
C       WMOID_LIST      WMO id's for the records in the read files.
C
C       NREC_LIST       Number of records (stations) in each read file.
C
C       ERROR_FLAG      Error handling flag.
C
        integer max_profilers
        parameter (max_profilers = 200)
 
        include 'prof_cdf_common.inc'
 
C
        character*21 error_text
        data error_text/'NetCDF file not open.'/
        integer i,j
C
        do i = 1, 4
C
C       Find the user's CDFID in our internal list.
C
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
C
              status = 0
              n_stations = nrec_list(i)
              if(n_stations.ne.0)then
C
                 do j = 1, nrec_list(i)
                    staname(j) = staname_list(j,i)
                    wmoid(j) = wmoid_list(j,i)
                 enddo
C
              endif
C
              go to 900
C
           endif                                ! End of found open file.
C
        enddo
C
        status = -5                             ! File not found.
C
900     if(status.ne.0.and.error_flag.ne.0)then
C
C       Use the non-default error handling for PROF_CDF errors.
C
                write(*,*)'PROF_CDF_GET_STATIONS:  '//error_text
                if(error_flag.eq.2)stop
        endif
C
        return
        end


        subroutine prof_i4_avg_wdw(i4_avg_wdw_sec, cdfid, istatus)
!  modified to pass in cdfid, and actually read file for data LW 8-27-98

        integer cdfid
	character*20 timestr
        character*13 attname
	integer    lenstr, sp_loc, units_loc, to_seconds, itime
        integer      tlen, ttype      ! attribute type and length

        include 'netcdf.inc'

C       netCDF file is opened, and accessed via cdfid
C	read "avgTimePeriod" global attribute into timestr 

        lenstr = len(timestr)  
        attname = 'avgTimePeriod'

C determine type of attribute 'avgTimePeriod'
        istatus = NF_INQ_ATTTYPE(cdfid,NF_GLOBAL,'avgTimePeriod',ttype)
        if (istatus .ne. NF_NOERR) then  !error retrieving info about avgTimePeriod
          istatus = 0
          return
        endif

        if (ttype .ne. NCCHAR) then   !read avgTimePeriod as an integer

          istatus = NF_GET_ATT_INT(cdfid,NF_GLOBAL,'avgTimePeriod',
     $                             itime)
	  if (istatus .ne. NF_NOERR) then   !error retrieving avgTimePeriod from file
            istatus = 0
            return
          endif
	  
C  units assumed to be minutes
	  to_seconds = 60

        else  !looking for a character string

C	  format of avgTimeString should be a number followed by a space and
C           then a lower case string of units, ie. "60 minutes" or "6 minutes"
C	    Verify that the units are minutes, seconds or hour, and convert
C	    the number string into an integer.  Return value is via i4_avg_wdw_sec
C  	    and is in seconds 

          istatus = NF_GET_ATT_TEXT(cdfid,NF_GLOBAL,'avgTimePeriod',
     $                              timestr)
	  if (istatus .ne. NF_NOERR) then   !error retrieving avgTimePeriod from file
            istatus = 0
            return
          endif

	  sp_loc = index(timestr,' ')  !everything to the left should be number
	
C	  convert numerical string to number
101	  format(i1)
102	  format(i2)
103	  format(i3)
104	  format(i4)
          itime = -1 ! check at end to see if set

	  if (sp_loc .eq. 2) then !one digit number
 	    read(timestr(1:1),101)itime
          endif

	  if (sp_loc .eq. 3) then !two digit number
 	    read(timestr(1:2),102)itime
          endif

	  if (sp_loc .eq. 4) then !three digit number
 	    read(timestr(1:3),103)itime
          endif

	  if (sp_loc .eq. 5) then !four digit number
 	    read(timestr(1:4),104)itime
          endif

	  if (itime .eq. -1) then   ! couldn't convert timestr to itime, return error
            write(6,*)
     1       ' Error: couldnt convert avgTimePeriod string to integer: '  
     1       , timestr     
	    istatus = 0
            return
	  endif

C	  determine units in timestr
          to_seconds = 0 !will check after looking for units to see if set
   	  units_loc = index(timestr,'minute')
          if (units_loc .gt. 0) to_seconds = 60
 
	  units_loc = index(timestr,'second')
          if (units_loc .gt. 0) to_seconds = 1
 
	  units_loc = index(timestr,'hour')
          if (units_loc .gt. 0) to_seconds = 3600

        endif

	if (to_seconds .eq. 0) then   ! couldn't identify units, return error
          write(6,*)
     1      ' Error: couldnt decode avgTimePeriod from file ', timestr       
	  istatus = 0
          return
        else
          i4_avg_wdw_sec = itime * to_seconds
          istatus = 1
	endif
        
        return
        end
