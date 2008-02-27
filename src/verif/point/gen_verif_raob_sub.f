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

      subroutine gen_verif_raob_sub(model_dir,i4time_sys
     1                             ,nl_dir, ni, nj, nk, stdLON
     1                             ,n,balance, r_missing_data
     1                             ,max_ht_m_proc, max_verif
     1                             ,istatus)

      implicit none

      character*9	a9_time_laps    !Input:  yyjjjhhmm from systime.dat
      character*9       a9_time_raob    !adjust time corresponding to i4time_raob
      character*256     model_dir	!Input:  location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time_sys      !Input:  i4time of LAPS analysis to verify (from systime.dat)
      integer         i4time_raob     !i4time of raob to process
      integer         i4time_laps     !internal i4time = to the analysis time
      integer           max_verif       !Input:  maximum # of verification types
      integer           ni, nj, nk	!Inputs: i, j and k grid dimensions
      real 		stdLON		!Input:  standard Longitude
      integer		balance         !Input:  1= processing balance; 0= not processing balance
      real		r_missing_data	!Input:  value used from LAPS for missing data
      real            max_ht_m_proc   !maximum height(m) to process up to
      integer           n               !Input:  reference to which verification type
      integer		istatus		!return value from subroutine

      character*9	a9time_raob(2)
      character*25      fileNames(1)
      integer		num_ret, i4timeraob(2)
      integer         i4time_raob_latest, i4time_raob_earliest
      integer		laps_cycle_time, i4_raob_window
      integer           min
      character*2	c_hr
      character*256	nl_dir		!directory where verify_raob.nl located
      integer		i4time_prev, n_raob_files
      integer           i4time_later
      character*256     raob_dir  	!path of raob file to read
      character*256     raob_fname	!filename of raob file to read
      character*256     output_dir  	!path for output file
      character*256     output_fname	!path and name of output file
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real            laps_levels_pa(nk) !laps pressure levels
      real            laps_levels_mb(nk) !laps pressure levels
      real            min_pres_mb_proc  !minimum pressure(mb) to process up to

      logical           lrunbal

      character*1       type_obs
      character*3       cpads_type
      character*8       c8_project
      character*150     path_to_raw_profiler
      character*150     path_to_raw_sounding
      integer  		raob_process_lag 
      integer  		raob_process_lag_bal
      integer           adv_anal_by_t_min
      integer           i_advanal_sec
      integer           n_verif
      integer           model_len
      character*150     verif_output_dir(max_verif)
c     character*150     verif_output_ext
c     character*150     verif_output_bal
      real		verif_missing_data
      logical           lexist
      integer		dir_len, ext_len, bal_len, i,ihr,lenn
C
C     BEGIN
C
C     read verify.nl to get raob directory and output directory
      write(6,*) 'in gen_sub'
      call read_verif_nl(type_obs,path_to_raw_profiler
     1                  ,path_to_raw_sounding, raob_process_lag
     1                  ,raob_process_lag_bal 
     1                  ,max_verif, verif_output_dir
     1                  ,verif_missing_data, n_verif, istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error in read_verif_nl '
        istatus = 0
        return
      endif

      i4time_laps=i4time_sys
      call make_fnam_lp(i4time_laps,a9_time_laps,istatus)

C     Process raob verif raob_process_lag seconds behind and
C     make output_fname
      call s_len(model_dir,model_len)
      call s_len(verif_output_dir(n),dir_len)
c     call s_len(verif_output_ext,ext_len)
c     call s_len(verif_output_bal,bal_len)

      call get_c8_project(c8_project,istatus)
      call upcase(c8_project,c8_project)
      call get_balance_nl(lrunbal,adv_anal_by_t_min,cpads_type,istatus)
      if(istatus.ne.0)then
         print*,'error getting balance namelist'
         stop
      endif

      i_advanal_sec=adv_anal_by_t_min*60

      if(balance .eq. 1)then
         i4time_raob=i4time_laps-raob_process_lag_bal  !for airdrop raob_process_lag= -1800
         call make_fnam_lp(i4time_raob, a9_time_raob, istatus)
         output_fname = model_dir(1:model_len)//'verif/'//
     1verif_output_dir(n)(1:dir_len)//'/'//a9_time_raob//'_raob'
         if(c8_project.eq.'AIRDROP'.and.balance.eq.1)then
            i4time_laps=i4time_laps+i_advanal_sec
            call make_fnam_lp(i4time_laps,a9_time_laps,istatus)
         endif
      elseif(balance.eq.0.or.balance.eq.2)then
         i4time_raob=i4time_laps-raob_process_lag  !for airdrop = raob_process_lag= 0
         call make_fnam_lp(i4time_raob, a9_time_raob, istatus)
         output_fname = model_dir(1:model_len)//'verif/'//
     1verif_output_dir(n)(1:dir_len)//'/'//a9_time_raob//'_raob'
c        if(c8_project.eq.'AIRDROP')then
c           i4time_laps=i4time_laps-i_advanal_sec
c           call make_fnam_lp(i4time_laps,a9_time_laps,istatus)
c        endif
      endif

      if(c8_project.ne.'AIRDROP')then
         i4time_laps=i4time_raob
         call make_fnam_lp(i4time_laps, a9_time_laps, istatus)
      endif

      write(6,*) 'Processing Analysis Time: ',a9_time_laps,
     1'  ',i4time_laps
  
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
      call get_pres_1d(0,nk,laps_levels_pa,istatus)
      if (istatus .ne. 1) then
        write(6,*)' Error getting LAPS pressure levels'
        istatus = 0
        return
      endif

C     convert pressure levels to mb
      do i = 1, nk
        laps_levels_mb(i) = laps_levels_pa(i)/100.
      enddo

C     set min_pres_mb_proc
      if (laps_levels_mb(1) .lt. laps_levels_mb(nk)) then
        min_pres_mb_proc = laps_levels_mb(1)
      else
        min_pres_mb_proc = laps_levels_mb(nk)
      endif

C     determine which raob file(s) to look for data in
C Caution: this likely will not work for cycle times other
C          than hourly or half hourly; but other times will
C          require some method to round to the nearest "top-of-hour".
      c_hr = a9_time_laps(6:7)
      read(a9_time_laps(8:9),'(i2)')min
      if(c_hr(1:1).eq.' ')c_hr(1:1)='0'
      a9time_raob(1) = a9_time_raob(1:5)//c_hr//'00'
      n_raob_files = 1
      call i4time_fname_lp(a9time_raob(1),i4timeraob(1), istatus)

      if(min.ge.30)then
         i4time_later=i4time_raob+min*60
         read(c_hr,'(i2)')ihr
         if(ihr.lt.23)then
            ihr=ihr+1
         else
            ihr=0
         endif
         write(c_hr,'(i2.2)')ihr
         call make_fnam_lp(i4time_later,a9time_raob(1),istatus)
         call i4time_fname_lp(a9time_raob(1),i4timeraob(1), istatus)
      endif

C     set i4_raob_window to laps_cycle_time
      call get_laps_cycle_time(laps_cycle_time,istatus)
      i4_raob_window = laps_cycle_time

C     Define limits of RAOB release times we are interested in
C     i4time_raob_latest =   i4timeraob(2) + i4_raob_window
C     i4time_raob_earliest = i4timeraob(1) - i4_raob_window
C     write(6,*) 'Raob time window= ',i4time_raob_earliest,
C    1           ' - ',i4time_raob_latest
 
C     set raob_dir
      raob_dir = path_to_raw_sounding
      call s_len(raob_dir,dir_len)

C     Loop through raob files

      do i = 1,n_raob_files

C       see if raob file is there
        raob_fname = a9time_raob(i)//'0300o'

c       call get_pres_1d(0,nk,laps_levels_pa,istatus)
c this call corrupts the pressures array in get_pres_1d.
c       call getFileNames_c(raob_dir,fileNames, num_ret,
c    1                      raob_fname,istatus)
c       call get_pres_1d(0,nk,laps_levels_pa,istatus)
c       if ((istatus.eq.1).and.(num_ret.eq.1)) then

        raob_fname = raob_dir(1:dir_len)//a9time_raob(i)//'0300o'
        call s_len(raob_fname,lenn)
        write(6,*) 'Looking for raob file ',raob_fname(1:lenn)
        inquire(file=raob_fname(1:lenn),exist=lexist)

        if(lexist)then

C     Define limits of RAOB release times we are interested in
          i4time_raob_latest =   i4timeraob(i) + i4_raob_window
          i4time_raob_earliest = i4timeraob(i) - i4_raob_window
          write(6,*) 'Raob time window= ',i4time_raob_earliest,
     1           ' - ',i4time_raob_latest


C         generate raob output file
          call get_raob_pairs(raob_fname,model_dir,i4time_laps
     1          ,i4time_raob,i4time_raob_earliest,i4time_raob_latest,
     1                        output_fname, nl_dir, ni, nj,
     1                        nk, lats, lons, stdLON, 
     1                        laps_levels_mb, laps_levels_pa,
     1                        max_ht_m_proc, min_pres_mb_proc,
     1                        balance, r_missing_data, 
     1                        verif_missing_data, istatus)

          if(istatus .ne. 1)then
            write(6,*)' Error in get_raob_pairs '
          endif
        else
          write(6,*) 'Raob file ',raob_fname,' not found.'
        endif
      enddo

      return
      end
!2.........................................................................
