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
      subroutine write_laps_lfm(reftime,valtime,datadir,cdldir, 
     1   ext,imax,jmax,
     1   kmax,kdim,var,lvl,lvl_coord,units,comment,n_levels,
     1   cdl_levels,data,istatus)

C**********************************************************************
C
C      This file contains the following FORTRAN subroutines:
C            write_laps_data
C
C      The write_laps_data subroutine reads the following FORTRAN
C      subroutines from the readlapsdata.f file:
C            cvt_fname_v3
C
C      The write_laps_data subroutine reads the following C subroutines
C      from the rwl_v3.c file:
C            write_cdf_v3
C
C**********************************************************************
C
C      Subroutine WRITE_LAPS_LVLS
C
C      Author:    John Snook
C      Modified:  Copied writelapsdata.f code and  6/98 Linda Wharton
C                   added parameter cdl_levels to
C                   make the program independent of
C                   nest7grid.parms and static.nest7grid
C
C      Writes data in arrays DATA and COMMENT to the netCDF file name
C      specified by I4TIME, DIR and EXT.  The data in VAR, LVL, LVL_COORD,
C      IMAX, JMAX, KMAX, KDIM and UNITS are stored into the netCDF file
C      when it is created.  ISTATUS is returned.
C
C**********************************************************************
C
        implicit  none
C
      integer      reftime,              !INPUT I4time of model run 
     1               valtime,              !INPUT I4time data is valid
     1               imax,jmax,kmax,       !INPUT # cols, # rows, # fields
     1               kdim,                 !INPUT K dimension of DATA array
     1               lvl(kdim),            !INPUT Level of each field 
     1               istatus               !OUTPUT

      real         data(imax,jmax,kdim)  !INPUT Raw data to be written
      character*(*)  datadir               !INPUT Directory to be written to
      character*(*)  cdldir                !INPUT Directory containing cdls
      character*(*)  ext                   !INPUT File name ext
      character*(*)  var(kdim)             !INPUT 3 letter ID of each field
      character*(*)  lvl_coord(kdim)       !INPUT Vertical coordinate of fields
      character*(*)  units(kdim)           !INPUT units of each field
      character*(*)  comment(kdim)         !INPUT Comments for each field
C
      integer      flag,                 !Print flag (1 = off)
     1               i_reftime,            !UNIX time of data
     1               i_valtime,            !UNIX time of data
     1               error(2),
     1               i,j,n7g_nx, n7g_ny,
     1               lgfc,
     1               fn_length,
     1               var_len,
     1               comm_len,
     1               ext_len,
     1               asc_len,
     1               lvl_coord_len,
     1               units_len,
     1               cdl_path_len,
     1               stat_len,
     1               max_levels,
     1               n_levels,
     1               called_from,          !0=FORTRAN, 1=C
     1               append                !0=no, 1=yes
C
      parameter (max_levels=100)
      real         pr(max_levels),       !pressures read from get_pres_1d
     1               cdl_levels(n_levels)
C
      character*5    fcst_hh_mm
      character*9    gtime
      character*150 file_name
      character*150  static_path
      character*9    laps_dom_file
      character*24   asctime
      character*40   v_g
      character*256  syscmd
      logical outfile_defined
C
      common         /prt/flag
C
      include 'grid_fname.cmn'
C
C-------------------------------------------------------------------------------
C
      error(1)=1
      error(2)=0
C 
C ****  Check if ext is fua or fsf or pbl
      if ((ext .eq. 'fua') .or. (ext .eq. 'fsf') .or.
     1 (ext .eq. 'pbl'))  then
C ****  Skip get_config 
      else
C
C ****  call get_config to read nest7grid.parms
C
        call get_config(istatus)

        call s_len(grid_fnam_common,lgfc)
        if (lgfc .gt. 0) then
          laps_dom_file = grid_fnam_common(1:lgfc)
        else
          laps_dom_file = grid_fnam_common
          write(6,*) 'Domain name not retrieved by get_config'
          write(6,*) 'Navigation info will not be written to',
     1' output file'
        endif

        call get_laps_dimensions(n_levels,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical domain dimension'
           return
        endif

        call get_grid_dim_xy(n7g_nx, n7g_ny, istatus)
        if (istatus .ne. 1) then
            write(6,*) 'return get_grid_dim_xy, status: ', istatus
            return
        endif

        call get_vertical_grid(v_g,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical grid'
           return
        endif

	if (v_g .ne. 'PRESSURE') goto 920

        call get_pres_1d(valtime,n_levels,pr,istatus)
        do j = 1,n_levels
          pr(j)=pr(j)/100.
        enddo
C
C ****  Various checks on input data.
C
        if (kmax .gt. kdim) then
          if (flag .ne. 1)
     1write (6,*) 'Illegal K dimension in DATA array...write aborted.'
          istatus=error(2)
          return
        endif
C
        if (imax .ne. n7g_nx) then
          if (flag .ne. 1)
     1write (6,*) 
     1'imax passed in does not match nest7grid.parms...write aborted.'
          istatus=error(2)
          return
        endif
C
        if (jmax .ne. n7g_ny) then
          if (flag .ne. 1)
     1write (6,*) 
     1'jmax passed in does not match nest7grid.parms...write aborted.'
          istatus=error(2)
          return
        endif

      endif
C
C ****  Get cdl_path***NOW PASSED IN AS CDLDIR 

C      call get_directory('cdl',cdl_path, cdl_path_len)
       call s_len(cdldir,cdl_path_len)
C ****  Get static_path
C
C      call get_directory('static',static_path, stat_len)
C
C ****  Specify file name
C
      call make_fnam_lp(reftime,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'Error converting i4time to file name...write aborted.'
        istatus=error(2)
        return
      endif

      call make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)
C
C ****  Create ascii time variables.
C
      call cv_i4tim_asc_lp(valtime,asctime,istatus)

      call s_len(ext, ext_len)

      call cvt_fname_v3(datadir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930

      called_from = 0    !called from FORTRAN
      append = 0         ! only one analysis time allowed per file

      var_len = len(var(1))
      comm_len = len(comment(1))
      lvl_coord_len = len(lvl_coord(1))
      units_len = len(units(1))
      asc_len = len(asctime)
      i_reftime = reftime - 315619200
      i_valtime = valtime - 315619200

      if(.true.)then 
!         Call ncgen here since it may not work from the C code
          if(imax*jmax*kmax .gt. 50000000)then
            syscmd = 'ncgen -v 2 -x -o '//trim(file_name)
     1               //' '//trim(cdldir)//trim(ext)//'.cdl'
          else
            syscmd = 'ncgen -x -o '//trim(file_name)
     1               //' '//trim(cdldir)//trim(ext)//'.cdl'
          endif
          write(6,*)' FORTRAN syscmd = ',trim(syscmd)
          call system(trim(syscmd))
          inquire(FILE=trim(file_name),EXIST=outfile_defined)
          if(outfile_defined .eqv. .false.)then
              write(6,*)' ERROR: output file does not exist'
              stop
          else
              write(6,*)
     1       ' output file successfully created via FORTRAN system call'
          endif
      else
          write(6,*)' use C code to run ncgen'
      endif
C
C **** write out netCDF file
C
      write(6,*)' Writing cdf file: ',file_name(1:fn_length)
      call write_cdf_v3 (file_name,ext,var,comment,asctime,cdldir, 
     1                   static_path,laps_dom_file, lgfc, fn_length,
     1                   ext_len,var_len, 
     1                   comm_len, asc_len, cdl_path_len, stat_len,
     1                   i_reftime, i_valtime,imax, jmax, kmax, kdim, 
     1                   lvl, data, pr, n_levels, cdl_levels,
     1                   called_from, append, istatus)
C
      if (istatus .gt. 0) goto 980
      IF (istatus .eq. -2) goto 940
      IF (istatus .eq. -3) goto 950
      IF (istatus .eq. -4) goto 960
      IF (istatus .eq. -5) goto 970
      IF (istatus .eq. -6) goto 990
C
C ****  Return normally.
C
        ISTATUS=ERROR(1)
999     RETURN
C
C ****  Error trapping.
C
920     IF (FLAG .NE. 1) THEN
          write(6,*) 'write_laps_lvls ABORTED!'
          write(6,*) ' LAPS will currently only work on a PRESSURE'
     1,' vertical grid'
          write(6,*) ' Make sure VERTICAL_GRID is set to PRESSURE'
     1,' in nest7grid.parms'
        ENDIF
        ISTATUS=ERROR(2)
        GOTO 999

930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...write aborted.'
        istatus=error(2)
        goto 999
C
940     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file to be written to...write abort
     1ed.'
        ISTATUS=ERROR(2)
        GOTO 999
C
950     IF (FLAG .NE. 1)
     1    write (6,*) 'Error in imax,jmax,or n_levels..write aborted'
        ISTATUS=ERROR(2)
        GOTO 999
C
960     IF (FLAG .NE. 1)
     1    write (6,*) 'Error writing data to file...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
970     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'Error writing header info into file...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
980     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'Some grids not written....could not convert LAPS variables.'
        ISTATUS=ERROR(2)
        GOTO 999
C
990     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'File already exists for analysis time...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
        END

C##########################################################################

