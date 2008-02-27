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

      subroutine gen_verif_prof_sub(model_dir, i4time, 
     1                              nl_dir, ni, nj, nk, stdLON,
     1                              n,balance, max_verif,
     1                              r_missing_data,istatus)

      implicit none

      character*9	a9_time
      character*256     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of LAPS/model file to read
      character*256	nl_dir		!directory where verify_prof.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real 		stdLON		!standard Longitude
      integer		balance
      real		r_missing_data	!value used from LAPS for missing data
      integer           max_verif       !Input: defines the maximum number of verification types
      integer           n_verif
      integer		istatus		!return value from subroutine

      character*256     prof_fname	!path and name of profiler file to read
      character*256     output_dir	!path for output file
      character*256     output_fname	!path and name of output file
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real            laps_levels_pa(nk) !laps pressure levels
      real            laps_levels_mb(nk) !laps pressure levels

      character*1       type_obs
      character*150     path_to_raw_profiler
      character*150     path_to_raw_sounding
      integer		raob_process_lag
      integer		raob_process_lag_bal
      character*150     verif_output_dir(max_verif)
c     character*150     verif_output_ext
c     character*150     verif_output_bal
      real		verif_missing_data
      integer		dir_len, model_len, bal_len, i,n
C
C     read verify.nl to get prof directory and output directory
      write(6,*) 'in gen_sub'
      call read_verif_nl(type_obs,path_to_raw_profiler,
     1 path_to_raw_sounding, raob_process_lag,raob_process_lag_bal,
     1                   max_verif, verif_output_dir,
     1                   verif_missing_data, n_verif, istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error in read_verif_nl '
        istatus = 0
        return
      endif

C     set prof_fname: this may have to move once we modify i4time and/or a9_time
c     based upon which verification we are performing (ie., balance = 0, or 1, or 2).
      call make_fnam_lp(i4time,a9_time,istatus)
      call s_len(path_to_raw_profiler,dir_len)
      prof_fname = path_to_raw_profiler(1:dir_len)//a9_time//'0100o'

C     make output_fname
      call s_len(model_dir,model_len)
      call s_len(verif_output_dir(n),dir_len)

c     call s_len(verif_output_ext,ext_len)
c     call s_len(verif_output_bal,bal_len)

c     if (balance .eq. 0) then

c        call s_len(verif_output_dir(1),dir_len)
c        output_fname = model_dir(1:model_len)//'verif/'//
c    1verif_output_dir(1)(1:dir_len)//'/'//a9_time//'_prof'
c        write(6,*) verif_output_dir(1)(1:dir_len)
c    1, '<',verif_output_dir(1)(1:dir_len),'>'
c     else
c-------------------------------------------------------------
         call s_len(verif_output_dir(n),dir_len)
         output_fname = model_dir(1:model_len)//'verif/'//
     1verif_output_dir(n)(1:dir_len)//'/'//a9_time//'_prof'
         write(6,*) verif_output_dir(2)(1:dir_len)
     1, '<',verif_output_dir(2)(1:dir_len),'>'
c-------------------------------------------------------------
c     endif
c replace the above with similar structures as in raob verif code.
c
c     call get_c8_project(c8_project,istatus)
c     call upcase(c8_project,c8_project)
c     call get_balance_nl(lrunbal,adv_anal_by_t_min,istatus)
c     if(istatus.ne.0)then
c        print*,'error getting balance namelist'
c        stop
c     endif

c     i_advanal_sec=adv_anal_by_t_min*60

c     if(balance .eq. 1)then
c        i4time_raob=i4time_laps-raob_process_lag_bal  !for airdrop raob_process_lag= -1800
c        call make_fnam_lp(i4time_raob, a9_time_raob, istatus)
c        output_fname = model_dir(1:model_len)//'verif/'//
c    1verif_output_dir(n)(1:dir_len)//'/'//a9_time_raob//'_raob'
c        if(c8_project.eq.'AIRDROP'.and.balance.eq.1)then
c           i4time_laps=i4time_laps+i_advanal_sec
c           call make_fnam_lp(i4time_laps,a9_time_laps,istatus)
c        endif
c     elseif(balance.eq.0.or.balance.eq.2)then
c        i4time_raob=i4time_laps-raob_process_lag  !for airdrop = raob_process_lag= 0
c        call make_fnam_lp(i4time_raob, a9_time_raob, istatus)
c        output_fname = model_dir(1:model_len)//'verif/'//
c    1verif_output_dir(n)(1:dir_len)//'/'//a9_time_raob//'_raob'
cc        if(c8_project.eq.'AIRDROP')then
cc           i4time_laps=i4time_laps-i_advanal_sec
cc           call make_fnam_lp(i4time_laps,a9_time_laps,istatus)
cc        endif
c     endif
c     if(c8_project.ne.'AIRDROP')then
c        i4time_laps=i4time_raob
c        call make_fnam_lp(i4time_laps, a9_time_laps, istatus)
c     endif

      write(6,*) 'Processing Analysis Time: ',a9_time    !_laps,
     1,'  ',i4time                                       !_laps

c
c careful if this is AIRDROP because in this case we want
c to verify analysis at t, balance at t+payload and backgrounds
c at either or both times (for comparison sake)

C     Read lat and lon data
      call read_static_grid(ni,nj,'LAT',lats,istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error getting LAPS LAT'
        istatus = 0
        return
      endif

      call read_static_grid(ni,nj,'LON',lons,istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error getting LAPS LON'
        istatus = 0
        return
      endif

C     Get LAPS/model pressure levels
      call get_pres_1d(i4time,nk,laps_levels_pa,istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error getting LAPS pressure levels'
        istatus = 0
        return
      endif

C     convert pressure levels to mb
      laps_levels_mb = laps_levels_pa/100.

C     generate profiler output file
      call get_prof_pairs(prof_fname, model_dir, i4time, 
     1                    output_fname, nl_dir, ni, nj,
     1                    nk, lats, lons, stdLON, 
     1                    laps_levels_mb, laps_levels_pa,
     1                    balance, r_missing_data, 
     1                    verif_missing_data, istatus)

      if(istatus .ne. 1)then
          write(6,*)' Error in get_prof_pairs '
      endif

      return
      end
!2.........................................................................
