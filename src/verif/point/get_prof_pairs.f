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

      subroutine get_prof_pairs(prof_fname, model_dir, i4time, 
     1                          output_fname, nl_dir, ni, nj,
     1                          nk, lats, lons, stdLON, 
     1                          laps_levels_mb,laps_levels_pa,
     1                          balance, r_missing_data, 
     1                          verif_missing_data, istatus)

      implicit none

      character*(*)     prof_fname	!path and name of profiler file to read
      character*(*)     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of LAPS/model file to read
      character*(*)     output_fname	!path and name of output file
      character*(*)	nl_dir		!directory where verify_prof.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real 		stdLON		!standard Longitude
      real            laps_levels_mb(nk) !laps pressure levels
      real            laps_levels_pa(nk) !laps pressure levels
      integer		balance
      real		r_missing_data
      real		verif_missing_data
      integer		istatus		!return value from subroutine

      integer		MAX_PROFS
      parameter         (MAX_PROFS=50)
      integer           numHTS
      parameter         (numHTS=64)
      integer		MAX_VAR
      parameter         (MAX_VAR=150)

      integer 		nf_status, ncid, varid
      integer           all_profs	! 1=use all profilers in domain
      integer           use_prof(MAX_PROFS), n_profs_use
      integer		kmax, kdim
      integer           lun, i, j, lend
      character*256     filename
      logical           l_eof
      integer 		wmoNum_use(MAX_PROFS)
      character*6       staName_use(MAX_PROFS)

      integer           nProfs, nHts(MAX_PROFS)
      integer		wmoStaNum(MAX_PROFS) 
      character*6       staName(MAX_PROFS) 
      real            staLat(MAX_PROFS), staLon(MAX_PROFS), 
     1                  staElev(MAX_PROFS), 
     1                  ri(numHTS,MAX_PROFS),
     1                  rj(numHTS,MAX_PROFS), 
     1                  rk(numHTS,MAX_PROFS)
      character*9	a9_time(MAX_PROFS)
      real            htMSL(numHTS, MAX_PROFS), 
     1                  uP(numHTS,MAX_PROFS), 
     1                  vP(numHTS,MAX_PROFS),
     1                  wP(numHTS,MAX_PROFS), 
     1                  baseLevels(numHTS)
      real            tSP(MAX_PROFS), rhSP(MAX_PROFS),
     1                  uSP(MAX_PROFS), vSP(MAX_PROFS)

      real            uI(numHTS,MAX_PROFS), 
     1                  vI(numHTS,MAX_PROFS),
     1                  wI(numHTS,MAX_PROFS), 
     1                  tSI(MAX_PROFS), rhSI(MAX_PROFS),
     1                  uSI(MAX_PROFS), vSI(MAX_PROFS)

      real            htM(ni,nj,nk) 
      real            tM(ni,nj,nk)
      real            uM(ni,nj,nk)
      real            vM(ni,nj,nk)
      real            omM(ni,nj,nk)
      real		tSM(ni,nj)
      real		rhSM(ni,nj)
      real		uSM(ni,nj)
      real		vSM(ni,nj)
C
C     BEGIN
C

      istatus = 1   !assume good return

C     Set up nHts(MAX_PROFS)
      do i = 1, MAX_PROFS
        nHts(i) = numHTS
      enddo

C     Read verify_prof.nl

      lun = 10
      call s_len(nl_dir,lend)
      filename = nl_dir(1:lend)//'/verif_prof.nl'
      open(lun,file=filename,status='old',err=900)

      l_eof = .false.
      n_profs_use = 0
      i = 1

50    format(i5,1x,a6)

      do while (.not.l_eof)
        read(lun,50,end=55,err=901)wmoNum_use(i), staName_use(i)
        i = i + 1
        goto 56
55      l_eof = .true.
56      continue
      enddo

      close(lun)
      n_profs_use = i - 1

      all_profs = 0	!set to 1 if use all profs in file

C     See if wmoNum_use(1) .eq. -1...if so, use all profilers in file
      if (wmoNum_use(1) .eq. -1) then
        all_profs = 1
      endif

C     Read profiler file "prof_fname" and return requested profiler data
      call get_prof_data(prof_fname, MAX_PROFS, numHTS, 
     1                   nProfs, a9_time,
     1                   stdLON, wmoStaNum, staName, 
     1                   staLat, staLon, staElev,
     1                   htMSL, tSP, rhSP, uSP,
     1                   vSP, uP, vP, wP, 
     1                   verif_missing_data, istatus)
       
      if (istatus .ne. 1) then
        write(6,*) ' Unable to read profiler file: ',prof_fname
        istatus = 0
        return
      endif

C     Set which profilers to pull LAPS/model data from
      if (all_profs .ne. 1) then   !set use_prof(i) to 1 if profiler is on list
        do i = 1, n_profs_use
          do j = 1, nProfs
            if ((wmoNum_use(i) .eq. wmoStaNum(j)) .or.
     1           (staName_use(i) .eq. staName(j))) then
              use_prof(j) = 1    !set to 1 if current profiler is on list to use
            else
              use_prof(j) = 0  
            endif
          enddo
        enddo
      else
        do i = 1, nProfs
          use_prof(i) = 1    !set to 1 if current profiler is on list to use
        enddo
      endif

C     get LAPS/model data
      call get_model_data(i4time, model_dir, ni, nj, nk, balance, 
     1                    laps_levels_mb, htM, uM, vM, omM, tM, 
     1                    tSM, rhSM, uSM, vSM, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' Unable to read model data '
        istatus = 0
        return
      endif

c     write(6,*) 'uSM(121,1)=',uSM(121,1)
c     write(6,*) 'vSM(121,1)=',vSM(121,1)
c     write(6,*) 'rhSM(121,1)=',rhSM(121,1)
c     write(6,*) 'tSM(121,1)=',tSM(121,1)
c     write(6,*) 
c     write(6,*) 'uSM(70,44)=',uSM(70,44)
c     write(6,*) 'vSM(70,44)=',vSM(70,44)
c     write(6,*) 'rhSM(70,44)=',rhSM(70,44)
c     write(6,*) 'tSM(70,44)=',tSM(70,44)
c     write(6,*) 
c     write(6,*) 'uSM(41,7)=',uSM(41,7)
c     write(6,*) 'vSM(41,7)=',vSM(41,7)
c     write(6,*) 'rhSM(41,7)=',rhSM(41,7)
c     write(6,*) 'tSM(41,7)=',tSM(41,7)
c     write(6,*) 

C     Loop through use_prof and pull model data if use_prof(i) .eq. 1
      do i = 1, nProfs
        if (use_prof(i) .eq. 1) then

C         Determine ri, rj, rk for ob in the model domain
          call get_rijk_ob_ht(i, MAX_PROFS, numHTS, 
     1                     laps_levels_pa, nHts,
     1                     staLat(i), staLon(i), lats,
     1                     lons, ni, nj, nk, htMSL, htM, 
     1                     ri, rj, rk, istatus)

          if (istatus .ne. 1) then
            use_prof(i) = 0
            if (istatus .eq. -1) then
            else
              write(6,*) ' Unable to determine ri,rj,rk for profiler ',
     1                 wmoStaNum(i), ' ',staName(i)
            endif
          else
C           Pull values interpolated from grid to ob
            call interp_to_ob(i, MAX_PROFS, numHTS,nHts,
     1                        ni, nj, nk, ri, rj, rk, uM, 
     1                        vM, omM, tM, tSM, rhSM, uSM,
     1                        vSM, uI, vI, wI, tSI, rhSI,
     1                        uSI, vSI, r_missing_data,
     1                        verif_missing_data, istatus) 

            if (istatus .ne. 1) then
              use_prof(i) = 0
              write(6,*) ' Unable to interp obs for profiler ',
     1                   wmoStaNum(i), ' ',staName(i)
            endif

          endif
        endif
      enddo

C     re-calc n_profs_use based on interpolation errors
      n_profs_use = 0
      do i = 1, nProfs
        if (use_prof(i) .eq. 1) n_profs_use = n_profs_use + 1
      enddo
       

C     write output_file

      call write_verif_prof(output_fname, MAX_PROFS, numHTS,
     1                      nHts, n_profs_use, use_prof, 
     1                      nProfs, wmoStaNum,
     1                      staName, staLat, staLon, staElev,
     1                      htMSL, ri, rj, rk, uP, vP, wP, tSP,
     1                      rhSP, uSP, vSP, uI, vI, wI, tSI,
     1                      rhSI, uSI, vSI, a9_time, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' Error writing profiler verif file: ',
     1              output_fname
      endif

      goto 999

900   print*,'error opening namelist file ', filename
      goto 999

901   print*,'error reading namelist file ', filename
      goto 999


999   continue

      return
      end
!1.........................................................................
      subroutine get_prof_data(infile, MAX_PROFS, numHTS, 
     1                         nProfs, a9_time,
     1                         stdLON, wmoStaNum, staName, 
     1                         staLat, staLon, staElev,
     1                         htMSL, tSfc, rhSfc, uSfc,
     1                         vSfc, u, v, w, 
     1                         verif_missing_data, istatus)

      implicit none
      include 'netcdf.inc'

      character*256	infile
      integer  		MAX_PROFS, numHTS, nProfs
      character*9       a9_time(MAX_PROFS)
      real		stdLON
      integer  		wmoStaNum(MAX_PROFS) 
      character*6       staName(MAX_PROFS) 
      real            staLat(MAX_PROFS), staLon(MAX_PROFS), 
     1                  staElev(MAX_PROFS), PI
      real            htMSL(numHTS, MAX_PROFS), 
     1                  u(numHTS,MAX_PROFS), v(numHTS,MAX_PROFS),
     1                  w(numHTS,MAX_PROFS), 
     1 			verif_missing_data,
     1                  baseLevels(numHTS)
      real            tSfc(MAX_PROFS), rhSfc(MAX_PROFS),
     1                  wsSfc(MAX_PROFS), wdSfc(MAX_PROFS),
     1                  uSfc(MAX_PROFS), vSfc(MAX_PROFS)
      real*8		timeObs(MAX_PROFS), prev_timeObs
      real		prof_missing

      integer           istatus, i, j, i4time
      integer           lenf
      character*9	prev_a9time
      integer           nf_fid, nf_vid, start_1(1), count_1(1),
     1                  start_2(2), count_2(2), nf_status

!..........................................................................
C     BEGIN

      istatus = 1   !assume good return

C  Open input file for reading
C
      call s_len(infile,lenf)
      print*,'Reading ',infile(1:lenf)
      nf_status = NF_OPEN(infile,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',infile
        istatus = 0
        return
      endif
C
C  Read number of profilers in file and number of levels
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nProfs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
        istatus = 0
        return
      endif

      if (nProfs .gt. MAX_PROFS) then
        print *,' MAX_PROFS too small to hold profilers in file.'
        istatus = 0
        return
      endif
C
C  Fill baseLevels
C
      start_1(1) = 1
      count_1(1) = 28

      nf_status = NF_INQ_VARID(nf_fid, 'levels', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid levels'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                            count_1, baseLevels(1))
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading levels 1:28'
        istatus = 0
        return
      endif

      start_1(1) = 37
      count_1(1) = 36

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                            count_1, baseLevels(29))
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading levels 37:72'
        istatus = 0
        return
      endif
C
C  Set start_1 and count_1 for reading 1D variables 
      start_1(1) = 1
      count_1(1) = nProfs
C
C  Read wmoStaNum 
C
      nf_status = NF_INQ_VARID(nf_fid, 'wmoStaNum', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid wmoStaNum'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_INT(nf_fid, nf_vid, start_1,
     1                            count_1, wmoStaNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading wmoStaNum'
        istatus = 0
        return
      endif
C
C  Read staName
C
      start_2(2) = 1
      count_2(2) = nProfs
      start_2(1) = 1
      count_2(1) = 6
C
      nf_status = NF_INQ_VARID(nf_fid, 'staName', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid staLat'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_TEXT(nf_fid, nf_vid, start_2, 
     1                             count_2, staName)

      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading staName'
        istatus = 0
        return
      endif

      do i = 1, nProfs
        staName(i)(6:6) = ' '
      enddo
C
C  Read staLat
C
      nf_status = NF_INQ_VARID(nf_fid, 'staLat', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid staLat'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, staLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading staLat'
        istatus = 0
        return
      endif
C
C  Read staLon
C
      nf_status = NF_INQ_VARID(nf_fid, 'staLon', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid staLon'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, staLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading staLon'
        istatus = 0
        return
      endif
C
C  Read staElev
C
      nf_status = NF_INQ_VARID(nf_fid, 'staElev', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid staElev'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, staElev)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading staElev'
        istatus = 0
        return
      endif
C
C  Read timeObs
C
      nf_status = NF_INQ_VARID(nf_fid, 'timeObs', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid timeObs'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_DOUBLE(nf_fid, nf_vid, start_1,
     1                             count_1, timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading timeObs'
        istatus = 0
        return
      endif

C     convert timeObs to a9_time string
      prev_timeObs = 0.
      prev_a9time = ' '
      do i = 1, nProfs
        if (timeObs(i) .eq. prev_timeObs) then
          a9_time(i) = prev_a9time
        else
          prev_timeObs = timeObs(i)
          i4time = nint(timeObs(i))  + 315619200
          call make_fnam_lp(i4time,prev_a9time,istatus)
          if(istatus .ne. 1)then
            prev_a9time = '         '
          endif
          a9_time(i) = prev_a9time
        endif
      enddo
C
C  Read surface temperature
C
      nf_status = NF_INQ_VARID(nf_fid, 'temperature', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid temperature'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, tSfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading temperature'
        istatus = 0
        return
      endif
C 
C  Read fill value for numerical variables
C
      nf_status = NF_GET_ATT_REAL(nf_fid, nf_vid, 
     1                            '_FillValue', 
     1                            prof_missing)
      
      do j = 1, nProfs
        if (tSfc(j) .eq. prof_missing) 
     1     tSfc(j) = verif_missing_data
      enddo
C
C  Read surface RH
C
      nf_status = NF_INQ_VARID(nf_fid, 'relHumidity', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid relHumidity'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, rhSfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading relHumidity'
        istatus = 0
        return
      else
        do j = 1, nProfs
          if (rhSfc(j) .gt. 100.0)
     1       rhSfc(j) = verif_missing_data
        enddo
      endif
C
C  Read surface wind speed
C
      nf_status = NF_INQ_VARID(nf_fid, 'windSpeedSfc', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid windSpeedSfc'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, wsSfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading windSpeedSfc'
        istatus = 0
        return
      else
        do j = 1, nProfs
          if (wsSfc(j) .eq. prof_missing) 
     1       wsSfc(j) = verif_missing_data
        enddo
      endif
C
C  Read surface wind direction 
C
      nf_status = NF_INQ_VARID(nf_fid, 'windDirSfc', nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'varid windDirSfc'
        istatus = 0
        return
      endif

      nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_1,
     1                             count_1, wdSfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading windDirSfc'
        istatus = 0
        return
      else
        do j = 1, nProfs
          if (wdSfc(j) .eq. prof_missing) 
     1       wdSfc(j) = verif_missing_data
        enddo
      endif
C
C  Loop through nProfs to read 1:28 and 37:72 of u, v, w
C
      do i = 1, nProfs
        start_2(2) = i
        count_2(2) = 1
        start_2(1) = 1
        count_2(1) = 28
C
C  Read u values
C
        nf_status = NF_INQ_VARID(nf_fid, 'uComponent', nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'varid uComponent'
          istatus = 0
          return
        endif

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, u(1,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading u 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (u(j,i) .eq. prof_missing) then
            u(j,i) = verif_missing_data
          else
            if (abs(u(j,i)) .gt. 200) then
              u(j,i) = verif_missing_data
              write(6,*) '*** ',nProfs, i, j, u(j,i) 
            endif
          endif
        enddo
        
        start_2(1) = 37
        count_2(1) = 36

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, u(29,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading u 37:72'
          istatus = 0
          return
        endif
 
c       if ((i .eq. 3) .or. (i .eq. nProfs)) then
c         write(6,*) 'u bad? 54,3:',u(54,3),' 55,3:',u(55,3)
c       endif

        do j = 29, 64
          if (u(j,i) .eq. prof_missing) then
            u(j,i) = verif_missing_data
          else
            if (abs(u(j,i)) .gt. 200) then
              write(6,*) '*** ',nProfs, i, j, u(j,i) 
              u(j,i) = verif_missing_data
            endif
          endif
        enddo
C
C  Read v values
C
        nf_status = NF_INQ_VARID(nf_fid, 'vComponent', nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'varid vComponent'
          istatus = 0
          return
        endif

        start_2(1) = 1
        count_2(1) = 28

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, v(1,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading v 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (v(j,i) .eq. prof_missing) then
            v(j,i) = verif_missing_data
          else
            if (abs(v(j,i)) .gt. 200) then
              write(6,*) '*** ',nProfs, i, j, v(j,i) 
              v(j,i) = verif_missing_data
            endif
          endif
        enddo

        start_2(1) = 37
        count_2(1) = 36

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, v(29,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading v 37:72'
          istatus = 0
          return
        endif

c       if ((i .eq. 3) .or. (i .eq. nProfs)) then
c         write(6,*) 'v bad? 54,3:',v(54,3),' 55,3:',v(55,3)
c       endif

        do j = 29, 64
          if (v(j,i) .eq. prof_missing) then
            v(j,i) = verif_missing_data
          else
            if (abs(v(j,i)) .gt. 200) then
              write(6,*) '*** ',nProfs, i, j, v(j,i) 
              v(j,i) = verif_missing_data
            endif
          endif
        enddo
C
C  Read w values
C
        nf_status = NF_INQ_VARID(nf_fid, 'wComponent', nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'varid wComponent'
          istatus = 0
          return
        endif

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, w(1,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading w 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (w(j,i) .eq. prof_missing) then
            w(j,i) = verif_missing_data
          else
            if (abs(w(j,i)) .gt. 20000.) then
              write(6,*) '*** ',nProfs, i, j, w(j,i) 
              w(j,i) = verif_missing_data
            endif
          endif
        enddo

        start_2(1) = 37
        count_2(1) = 36

        nf_status = NF_GET_VARA_REAL(nf_fid, nf_vid, start_2,
     1                               count_2, w(29,i))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading w 37:72'
          istatus = 0
          return
        endif

c     if ((i .eq. 15) .or. (i .eq. nProfs)) then
c       write(6,*) 'w bad? 29,15:',v(29,15)
c       write(6,*) 'w bad? 39,15:',v(39,15)
c       write(6,*) 'w bad? 45,15:',v(45,15)
c       write(6,*) 'w bad? 47,15:',v(47,15)
c       write(6,*) 'w bad? 49,15:',v(49,15)
c       write(6,*) 'w bad? 55,15:',v(55,15)
c       write(6,*) 'w bad? 57,15:',v(57,15)
c     endif

        do j = 29, 64
          if (w(j,i) .eq. prof_missing) then
            w(j,i) = verif_missing_data
          else
            if (abs(w(j,i)) .gt. 20000.) then
              write(6,*) '*** ',nProfs, i, j, w(j,i) 
              w(j,i) = verif_missing_data
            endif
          endif
        enddo

      enddo
C
C  Close input file
C
      nf_status = NF_CLOSE(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_CLOSE ',infile
        istatus = 0
        return
      endif
C
C  Fill htMSL for each profiler
C
      do i = 1, nProfs
        do j = 1, numHTS
          htMSL(j,i) = staElev(i) + baseLevels(j)
        enddo
      enddo
C
C  Convert wsSfc and wdSfc to uSfc and vSfc
C
      do i = 1, nProfs
        if ((wdSfc(i) .eq. verif_missing_data) .or.
     1      (wdSfc(i) .eq. verif_missing_data)) then
          uSfc(i) = verif_missing_data
          vSfc(i) = verif_missing_data
        else
          call disptrue_to_uvgrid(wdSfc(i), 
     1                            wsSfc(i),
     1                            uSfc(i),
     1                            vSfc(i),
     1                            stdLON)
          if (uSfc(i) .gt. 100) then
            uSfc(i) = verif_missing_data
            write(6,*) '*** ',nProfs, i, uSfc(i) 
          endif
          if (vSfc(i) .gt. 100) then
            write(6,*) '*** ',nProfs, i, vSfc(i) 
            vSfc(i) = verif_missing_data
          endif
        endif
      enddo

c     if ((i .eq. 15) .or. (i .eq. nProfs)) then
c       write(6,*) 'tSfc bad? 15:',tSfc(15)
c       write(6,*) 'rhSfc bad? 15:',rhSfc(15)
c       write(6,*) 'uSfc bad? 15:',uSfc(15)
c       write(6,*) 'vSfc bad? 15:',vSfc(15)
c     endif

      return
      end
!2.........................................................................
      subroutine get_model_data(i4time, model_dir, ni, nj, nk, 
     1                          balance, laps_levels,
     1                          ht, u, v, om, t, 
     1                          tS, rhS, uS, vS, istatus)

      implicit none

      integer		MAX_VAR
      parameter         (MAX_VAR=150)

      integer		i4time		!i4time of LAPS/model file to read
      character*(*)     model_dir	!location of model data directories
      integer           ni, nj, nk	!i, j and k grid dimensions
      integer		balance
      real            ht(ni,nj,nk)     !data to be read
      real            t(ni,nj,nk)     !data to be read
      real            u(ni,nj,nk)     !data to be read
      real            v(ni,nj,nk)     !data to be read
      real            om(ni,nj,nk)     !data to be read
      real		tS(ni,nj)
      real		rhS(ni,nj)
      real		uS(ni,nj)
      real		vS(ni,nj)
      real		pS(ni,nj)
      real            make_rh
      integer		istatus		!return value from subroutine

      character*256     dir_in
      real            laps_levels(nk) !laps pressure levels
      character*3       var(MAX_VAR)
      character*4       lvl_coord(MAX_VAR)      !Vertical coordinate for each field
      character*10      units(MAX_VAR)    	!units of each field
      character*125     comment(MAX_VAR)  	!Comments for each field
      integer           lvl(nk)
      integer           lend, i,j, nkS
      integer           dir_len
      character*31      ext
      character*3       csubdir
      integer           lt
      integer           i4time_fcst
      integer           i4time_init
      character*20      cmdltype
      character*9       a9_time
C
C     BEGIN
C
      istatus = 1   !assume good return

      call make_fnam_lp(i4time,a9_time,istatus)
      call bkgd_wgi(a9_time,i4time_init,i4time_fcst
     +,csubdir,cmdltype,balance,istatus)
      if(istatus.ne.1)then
         print*,'Failure in bkgd_wgi to get model bkgd time'
         return
      endif

      call s_len(model_dir, dir_len)
      lend=dir_len
      dir_in=model_dir
      if(dir_in(lend:lend).ne.'/')then
         lend=lend+1
         dir_in(lend:lend)='/'
      endif

      if(balance .eq. 1)then
        dir_in = model_dir(1:lend)//'balance/'
      elseif(balance .eq. 2)then      !this is the background;
c                                      must be that used in the analysis.
          if(csubdir.eq.'fua')then
             call s_len(cmdltype,lt)
      dir_in=model_dir(1:lend)//csubdir//'/'//cmdltype(1:lt)//'/'
          else
             dir_in=model_dir(1:lend)//'lga/'
          endif

      endif

C     Read LAPS/model grid  ht
      ext = 'lt1'
      var = 'HT'
      lvl = laps_levels

      if(balance .ne. 2)then

         dir_in=dir_in(1:lend)//'lt1/'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,ht,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading HT variable '
             return
         endif
         var = 'T3'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,t,istatus)

         if (istatus .ne. 1) then
             print*,' Error reading T3 variable '
             return
         endif

      else

         ext = csubdir
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,ht,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading HT variable '
             return
         endif
         var = 'T3'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,t,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading T3 variable '
             return
         endif

      endif

C     Read LAPS/model grid  u
      ext = 'lw3'
      var = 'U3'
      if(balance.ne.2)then

         dir_in=dir_in(1:lend)//'lw3/'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,u,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading U3 variable '
             return
         endif
         var = 'V3'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,v,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading V3 variable '
             return
         endif
         var = 'OM'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,om,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading OM variable '
             return
         endif

      else

         ext = csubdir
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,u,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading U3 variable '
             return
         endif
         var = 'V3'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,v,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading V3 variable '
             return
         endif
         var = 'OM'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,om,istatus)
         if (istatus .ne. 1) then
             print*,' Error reading OM variable '
             return
         endif

      endif

      print*,'LAPS 3D: t, ht, u, v, om successfully read'

C     Read LAPS/model grid  SFC u, v, t, rh

      if(balance.eq.0.or.balance.eq.1)then
         ext = 'lsx'
         dir_in = model_dir(1:dir_len)//'/lsx/'
         nkS = 1
         lvl(1) = 0
         var(1) = 'U'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nkS,
     1                     nkS,var,lvl,lvl_coord,
     1                     units,comment,uS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: u '
             return
         endif

         var(1) = 'V'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nkS,
     1                     nkS,var,lvl,lvl_coord,
     1                     units,comment,vS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: v '
             return
         endif

         var(1) = 'T'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nkS,
     1                     nkS,var,lvl,lvl_coord,
     1                     units,comment,tS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: t '
             return
         endif

         var(1) = 'RH'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nkS,
     1                     nkS,var,lvl,lvl_coord,
     1                     units,comment,rhS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: rh '
             return
         endif

      elseif(balance.eq.2)then

         if(csubdir.eq.'fua')then
      dir_in=model_dir(1:dir_len)//'fsf/'//cmdltype(1:lt)//'/'
            ext = 'fsf'
         else
            dir_in = model_dir(1:dir_len)//'/lgb/'
            ext = 'lgb'
         endif
         nkS = 1
         lvl(1) = 0
         var(1) = 'USF'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,uS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: usf'
             return
         endif

         var(1) = 'VSF'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,vS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: vsf'
             return
         endif

         var(1) = 'TSF'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,tS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: tsf'
             return
         endif

         var(1) = 'PSF'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,pS,istatus)
         if (istatus .ne. 1) then
             print*,' Unable to read model Surface variable: psf '
             return
         endif

         if(ext.eq.'fsf')then
            var(1) = 'RH'
            call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,rhS,istatus)
            if (istatus .ne. 1) then
                print*,' Unable to read model Sfc var: rh '
             return
            endif

         else

            var(1) = 'RSF'
            call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nkS,nkS,var,lvl,lvl_coord,
     1                  units,comment,rhS,istatus)
            if (istatus .ne. 1) then
                print*,' Unable to read model Sfc var: rsf '
                return
            endif

c note: pS is in pa so no need to multiply result by 100 (rh 0-100%)
            do j=1,nj
            do i=1,ni
               rhS(i,j)=make_rh(pS(i,j),tS(i,j)-273.15,rhS(i,j)*1000.
     +,-132.)
            enddo
            enddo

         endif

      endif

      write(6,*) 'LAPS Sfc t, rh, u, v successfully read'

      return
      end

!3.........................................................................
      subroutine get_rijk_ob_ht(whichIndex, MAX_PROFS, numHTS,
     1                          laps_levels, nHts,
     1                          lat_ob, lon_ob, lats,
     1                          lons, ni, nj, nk, htMSL, 
     1                          ht, ri, rj, rk, istatus)

      implicit none

      integer		whichIndex
      integer		MAX_PROFS, numHTS
      real            laps_levels(nk) !laps pressure levels
      integer		nHts(MAX_PROFS)
      real		lat_ob, lon_ob
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      integer           ni, nj, nk	!i, j and k grid dimensions
      real            htMSL(numHTS,MAX_PROFS), 
     1                  ht(ni,nj,nk),
     1                  ri(numHTS,MAX_PROFS),
     1                  rj(numHTS,MAX_PROFS), 
     1                  rk(numHTS,MAX_PROFS)
      integer		istatus

      integer		i, j, int_ri, int_rj
      real		height_to_zcoord3

C
C     BEGIN
C
      istatus = 1   !assume good return

C     Get ri, rj of ob 
      call latlon_to_rlapsgrid(lat_ob,lon_ob,lats,lons,ni,
     1                         nj,ri(1,whichIndex),
     1                         rj(1,whichIndex), istatus)
 
      if (istatus .eq. 1) then
        write(6,*)
        write(6,*) ' Profiler in domain: ', whichIndex
      else
c       write(6,*) ' Profiler out of domain: ',whichIndex
        istatus = -1
        return
      endif

C     Fill rest of ri,rj for height column of ob
      do i = 2, numHTS
        ri(i,whichIndex) = ri(1,whichIndex)
        rj(i,whichIndex) = rj(1,whichIndex)
      enddo

C     get integral ri, rj for height_to_zcoord3  
      int_ri = nint(ri(1,whichIndex))
      int_rj = nint(rj(1,whichIndex))

C     Calc rk for each level
      do i = 1, numHTS
        rk(i,whichIndex) = height_to_zcoord3(htMSL(i,whichIndex),
     1                                       ht,laps_levels,
     1                                       ni,nj,nk,int_ri,
     1                                       int_rj,istatus)

        if (istatus .ne. 1) then
          if (i .gt. 1) then
            nHts(whichIndex) = i - 1
            istatus = 1
            return
          else
            write(6,*) ' Unable to get rk for profiler ',whichIndex,
     1               '  height index= ',i,'maxht=',
     1               ht(int_ri,int_rj,nk),'pht=',htMSL(i,whichIndex)
            istatus = 0
            return
          endif
        endif
      enddo

      return
      end
!4.........................................................................
      subroutine write_verif_prof(output_fname, MAX_PROFS, numHTS,
     1                            nHts, n_profs_use, use_prof, 
     1                            nProfs, wmoStaNum,
     1                            staName, staLat, staLon, staElev,
     1                            htMSL, ri, rj, rk, uP, vP, wP, tSP,
     1                            rhSP, uSP, vSP, uI, vI, wI, tSI,
     1                            rhSI, uSI, vSI, a9_time, istatus)

      implicit none

      character*(*)     output_fname	!path and name of output file
      integer		MAX_PROFS, numHTS, nHts(MAX_PROFS)
      integer           n_profs_use, use_prof(MAX_PROFS), nProfs
      integer		wmoStaNum(MAX_PROFS) 
      character*6       staName(MAX_PROFS) 
      real            staLat(MAX_PROFS), staLon(MAX_PROFS), 
     1                  staElev(MAX_PROFS), 
     1                  htMSL(numHTS, MAX_PROFS), 
     1                  ri(numHTS,MAX_PROFS),
     1                  rj(numHTS,MAX_PROFS), 
     1                  rk(numHTS,MAX_PROFS),
     1                  uP(numHTS,MAX_PROFS), 
     1                  vP(numHTS,MAX_PROFS),
     1                  wP(numHTS,MAX_PROFS), 
     1                  tSP(MAX_PROFS), rhSP(MAX_PROFS),
     1                  uSP(MAX_PROFS), vSP(MAX_PROFS)
      real            uI(numHTS,MAX_PROFS), 
     1                  vI(numHTS,MAX_PROFS),
     1                  wI(numHTS,MAX_PROFS), 
     1                  tSI(MAX_PROFS), rhSI(MAX_PROFS),
     1                  uSI(MAX_PROFS), vSI(MAX_PROFS)
      character*9	a9_time(MAX_PROFS)
      integer		istatus

      integer		i, j

C
C     BEGIN
C
      istatus = 1   !assume good return

      write(6,*) output_fname

C     open output_fname
      open(1,file=output_fname,status='unknown',err=98)
      go to 99

98    write(6,*)' Error opening verif file: ',output_fname
      istatus = 0
      return

99    continue

C 100 writes: n_profs_use
100   format(i3)
C 101 writes: nHts+sfc,a9_time,wmoStaNum,staName,staLat,staLon,staElev
101   format(i3,1x,a9,1x,i7,1x,a6,1x,f7.3,1x,f8.3,1x,f7.0)
C 102 writes: staElev,ri,rj,uSP,uSI,vSP,vSI,tSP,tSI,rhSP,rhSI
102   format(f7.0,1x,2(f8.3,1x),9x,8(f6.1,1x))
C 103 writes: htMSL,ri,rj,rk,uP,uI,vP,vI,wP,wI
103   format(f7.0,1x,3(f8.3,1x),6(f6.1,1x))
C104  write: Header for sfc variables
104   format(1x,'Stnhgt',4x,'ri',7x,'rj',16x,'uSP',4x,'uSI',
     +4x,'vSP',4x,'vSI',4x,'tSP',4x,'tSI',3x,'rhSP',3x,'rhSI')
C 1010writes: Header for UA variables
105   format(4x,'Hgt',4x,'ri',7x,'rj',7x,'rk',8x,'uP',5x,
     +'uI',5x,'vP',5x,'vI',5x,'wP',5x,'wI')
C     write number of profilers into file
      write(1,100) n_profs_use
      
C     write each profiler out where use_prof(i) .eq. 1
      do i = 1, nProfs
        if (use_prof(i) .eq. 1) then
C         write header info
          write(1,101) nHts(i)+1,a9_time(i),wmoStaNum(i),
     1                 staName(i)(1:6),
     1                 staLat(i),staLon(i),staElev(i) 
C         write surface data
          write(1,104)
          write(1,102) staElev(i),ri(1,i),rj(1,i),uSP(i),uSI(i),
     1                 vSP(i),vSI(i),tSP(i),tSI(i),rhSP(i),rhSI(i)

          write(1,105)
          do j = 1, nHts(i)
            write(1,103) htMSL(j,i),ri(j,i),rj(j,i),rk(j,i),uP(j,i),
     1                   uI(j,i),vP(j,i),vI(j,i),wP(j,i),wI(j,i)
          enddo
        endif
      enddo

      close(1)

      return
      end
!5.........................................................................
      subroutine interp_to_ob(whichIndex, MAX_PROFS, numHTS,
     1                        nHts,
     1                        ni, nj, nk, ri, rj, rk, uM, 
     1                        vM, omM, tM, tSM, rhSM, uSM,
     1                        vSM, uI, vI, wI, tSI, rhSI,
     1                        uSI, vSI, r_missing_data,
     1                        verif_missing_data, istatus) 

      implicit none

      integer		whichIndex, MAX_PROFS, numHTS
      integer           nHts(MAX_PROFS)
      integer           ni, nj, nk	!i, j and k grid dimensions
      real            ri(numHTS,MAX_PROFS),
     1                  rj(numHTS,MAX_PROFS), 
     1                  rk(numHTS,MAX_PROFS),
     1                  uI(numHTS,MAX_PROFS), 
     1                  vI(numHTS,MAX_PROFS),
     1                  wI(numHTS,MAX_PROFS), 
     1                  tSI(MAX_PROFS), rhSI(MAX_PROFS),
     1                  uSI(MAX_PROFS), vSI(MAX_PROFS),
     1                  tM(ni,nj,nk), uM(ni,nj,nk),
     1                  vM(ni,nj,nk), omM(ni,nj,nk),
     1      		tSM(ni,nj), rhSM(ni,nj),
     1      		uSM(ni,nj), vSM(ni,nj)
      real		r_missing_data, verif_missing_data
      integer		istatus

      integer           int_ri, int_rj, i, extrap
      real		om, temp, pres
      real		pressure_of_rlevel
C
C     BEGIN
C
      istatus = 1   !assume good return

c     write(6,*) 'ri=',ri(1,whichIndex),' rj=',rj(1,whichIndex)
c     int_ri = nint(ri(1,whichIndex))
c     int_rj = nint(rj(1,whichIndex))
c     write(6,*) 'int_ri=',int_ri,' int_rj=',int_rj
c     write(6,*)

C     determine if extrapolation is needed 
      if (((ri(1,whichIndex) .lt. 1.0) .and. 
     1     (ri(1,whichIndex) .gt. 0.5)) 
     1    .or.
     1    ((ri(1,whichIndex) .lt. ni+0.5) .and. 
     1     (ri(1,whichIndex) .gt. ni))
     1    .or.
     1    ((rj(1,whichIndex) .lt. 1.0) .and. 
     1     (ri(1,whichIndex) .gt. 0.5))
     1    .or.
     1    ((rj(1,whichIndex) .lt. nj+0.5) .and. 
     1     (ri(1,whichIndex) .gt. nj))) then
        extrap = 1
      else
        extrap = 0
      endif

C     interpolate surface variable: uSI from uSM
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichIndex),
     1                              rj(1,whichIndex),ni,nj,
     1                              uSM,uSI(whichIndex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichIndex),rj(1,whichIndex),
     1                     ni,nj,uSM,uSI(whichIndex))
      endif

      if (uSI(whichIndex) .eq. r_missing_data)
     1  uSI(whichIndex) = verif_missing_data

c     write(6,*) 'uSM(int_ri,int_rj)=',uSM(int_ri,int_rj)
c     write(6,*) 'uSI(',whichIndex,')=',uSI(whichIndex)
c     write(6,*)

C     interpolate surface variable: vSI from vSM
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichIndex),
     1                              rj(1,whichIndex),ni,nj,
     1                              vSM,vSI(whichIndex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichIndex),rj(1,whichIndex),
     1                     ni,nj,vSM,vSI(whichIndex))
      endif

      if (vSI(whichIndex) .eq. r_missing_data)
     1  vSI(whichIndex) = verif_missing_data

c     write(6,*) 'vSM(int_ri,int_rj)=',vSM(int_ri,int_rj)
c     write(6,*) 'vSI(',whichIndex,')=',vSI(whichIndex)
c     write(6,*)

C     interpolate surface variable: rhSI from rhSM
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichIndex),
     1                              rj(1,whichIndex),ni,nj,
     1                              rhSM,rhSI(whichIndex),
     1                              istatus)
      else

        call bilinear_laps(ri(1,whichIndex),rj(1,whichIndex),
     1                     ni,nj,rhSM,rhSI(whichIndex))
      endif

      if (rhSI(whichIndex) .eq. r_missing_data)
     1  rhSI(whichIndex) = verif_missing_data

c     write(6,*) 'rhSM(int_ri,int_rj)=',rhSM(int_ri,int_rj)
c     write(6,*) 'rhSI(',whichIndex,')=',rhSI(whichIndex)
c     write(6,*)

C     interpolate surface variable: tSI from tSM
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichIndex),
     1                              rj(1,whichIndex),ni,nj,
     1                              tSM,tSI(whichIndex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichIndex),rj(1,whichIndex),
     1                     ni,nj,tSM,tSI(whichIndex))
      endif

      if (tSI(whichIndex) .eq. r_missing_data)
     1  tSI(whichIndex) = verif_missing_data

c     write(6,*) 'tSM(int_ri,int_rj)=',tSM(int_ri,int_rj)
c     write(6,*) 'tSI(',whichIndex,')=',tSI(whichIndex)
c     write(6,*)

C     loop through levels for u,v,w
      do i = 1, nHts(whichIndex)

C       interpolate 3D variable: uI from uM
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichIndex),
     1                                 rj(i,whichIndex),
     1                                 rk(i,whichIndex),
     1                                 ni,nj,nk, uM,
     1                                 uI(i,whichIndex),
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichIndex),rj(i,whichIndex),
     1                        rk(i,whichIndex),ni,nj,nk, uM,
     1                        uI(i,whichIndex))
        endif

        if (uI(i,whichIndex) .eq. r_missing_data)
     1    uI(i,whichIndex) = verif_missing_data

C       interpolate 3D variable: vI from vM
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichIndex),
     1                                 rj(i,whichIndex),
     1                                 rk(i,whichIndex),
     1                                 ni,nj,nk, vM,
     1                                 vI(i,whichIndex),
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichIndex),rj(i,whichIndex),
     1                        rk(i,whichIndex),ni,nj,nk, vM,
     1                        vI(i,whichIndex))
        endif

        if (vI(i,whichIndex) .eq. r_missing_data)
     1    vI(i,whichIndex) = verif_missing_data

C       interpolate 3D variable: om from omM
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichIndex),
     1                                 rj(i,whichIndex),
     1                                 rk(i,whichIndex),
     1                                 ni,nj,nk, omM,
     1                                 om,
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichIndex),rj(i,whichIndex),
     1                        rk(i,whichIndex),ni,nj,nk, omM,
     1                        om)
        endif

C       interpolate 3D variable: temp from tM
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichIndex),
     1                                 rj(i,whichIndex),
     1                                 rk(i,whichIndex),
     1                                 ni,nj,nk, tM,
     1                                 temp,
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichIndex),rj(i,whichIndex),
     1                        rk(i,whichIndex),ni,nj,nk, tM,
     1                        temp)
        endif

        if ((om .eq. r_missing_data) .or.
     1      (temp .eq. r_missing_data))  then
          wI(i,whichIndex) = verif_missing_data
        else
C         get pressure of rk
          pres = pressure_of_rlevel(rk(i,whichIndex))

C         calclulate wI(i,whichIndex) from om   w = -om*(R*temp/g*pres)
          wI(i,whichIndex) = -1*om*((287.04*temp)/(9.80815*pres))
        endif

      enddo

      return
      end
!6.........................................................................

