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

      program gen_verif_raob

      implicit none

      include 'grid_fname.cmn'

      integer           max_verif
      parameter        (max_verif=4)

      character*256     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of LAPS/model file to read
      character*9	a9_time
      character*150     verif_output_dir(max_verif) ! -->
                                        !type of verif  desired ... builds output directory
                                        !names (noBal, Bal, Bkgd, grid)
      character*256	nl_dir		!directory where verify_raob.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real 		stdLON		!standard Longitude
      integer		istatus		!return value from subroutines
      integer		len, balance,i
      real		r_missing_data
      character*2       cur_hr
      integer		i4_raob_window 
      real		max_ht_m_proc     !maximum height(m) to process up to

      character*1       type_obs
      character*150     path_to_raw_profiler
      character*150     path_to_raw_sounding
      integer           raob_process_lag_Bal
      integer           raob_process_lag
      integer           n_verif
      real              verif_missing_data
C
C     BEGIN
C

C     max height to process up to in meters
      max_ht_m_proc = 25000.0

C     do LAPS housekeeping

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
         write(6,*)' error in get_laps_config'
         goto 999
      endif

      call get_grid_dim_xy(ni,nj,istatus)

      if (istatus .ne. 1) then
         write (6,*) 'Error getting horizontal domain dimensions'
         go to 999
      endif

      call get_laps_dimensions(nk,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'Error getting vertical domain dimension'
         go to 999
      endif

C     get stdLON
      call get_standard_longitude(stdLON,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'Error getting standard longitude'
         go to 999
      endif

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1) then
        write(6,*) 'Error gettin r_missing_data value'
        goto 999
      endif

C     get current time from systime.dat
      call get_systime(i4time,a9_time,istatus)

      if(istatus .ne. 1) then
        write(6,*) 'Error from get_systime'
        go to 999
      endif

C DEBUG LW
c     i4time = 1346508000
c     a9_time = '022441400'

      call read_verif_nl(type_obs,path_to_raw_profiler,
     1 path_to_raw_sounding, raob_process_lag,raob_process_lag_bal,
     1                   max_verif, verif_output_dir,
     1                   verif_missing_data, n_verif, istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error in read_verif_nl '
        istatus = 0
        stop
      endif

      

C     get the directory for model_dir (lapsprd directory)
C     get the directory for nl_dir (static directory)

      call get_directory('lw3',model_dir,len)
      model_dir = model_dir(1:len-4)
      call get_directory(grid_fnam_common,nl_dir,len)

      do i=1,n_verif

         call downcase(verif_output_dir(i),verif_output_dir(i))
         if(verif_output_dir(i)(1:len).eq.'bal')then
            balance=1
         elseif(verif_output_dir(i)(1:len).eq.'nobal')then
            balance=0
         elseif(verif_output_dir(i)(1:len).eq.'bkgd')then
            balance=2
         endif
         

         call gen_verif_raob_sub(model_dir, i4time, 
     1                        nl_dir, ni, nj, nk, stdLON,
     1                        i, balance, r_missing_data, 
     1                        max_ht_m_proc, max_verif,
     1                        istatus)

         if(istatus .ne. 1) then
             write (6,*) 'Error in gen_verif_raob_sub'
         endif

      enddo

999   continue

      end
