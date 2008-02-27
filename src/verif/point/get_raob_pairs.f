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
      subroutine get_raob_pairs(raob_fname,model_dir,i4time_sys,i4time,
     1                          i4time_raob_earliest,i4time_raob_latest,
     1                          output_fname, nl_dir, ni, nj,
     1                          nk, lats, lons, stdLON, 
     1                          laps_levels_mb,laps_levels_pa,
     1                          max_ht_m_proc, min_pres_mb_proc,
     1                          balance, r_missing_data, 
     1                          verif_missing_data, istatus)

      implicit none

      character*(*)     raob_fname	!path and name of raob file to read
      character*(*)     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of LAPS/model file to read
      integer           i4time_sys
      integer         i4time_raob_latest, i4time_raob_earliest
      character*256     output_fname	!path and name of output file
      character*(*)	nl_dir		!directory where verify_raob.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real 		stdLON		!standard Longitude
      real            laps_levels_mb(nk) !laps pressure levels
      real            laps_levels_pa(nk) !laps pressure levels
      real            max_ht_m_proc   !maximum height(m) to process up to
      real            min_pres_mb_proc  !minimum pressure(mb) to process up to

      integer		balance
      real		r_missing_data, raob_missing_data
      real		verif_missing_data
      integer		istatus		!return value from subroutine

      integer		MAX_RAOBS
      parameter         (MAX_RAOBS=50)
      integer           maxM
      parameter         (maxM=22)
      integer           maxW
      parameter         (maxW=97)  !75+22, for mingle mand and sigW
      integer		maxT
      parameter         (maxT=177) !150+22 for mingle mand and sigT
      integer		MAX_HTS
      parameter		(MAX_HTS=275)

      integer 		nf_status, ncid, varid
      integer           all_raobs	! 1=use all raobs in domain
      integer           use_raob(MAX_RAOBS), n_raobs_use
      integer           lun, i, j, dir_len
      integer 	statusL(2,6)
      character*256     filename
      logical           l_eof

C     1d raob data
      integer		nRaobs, n_raobs_avail
      integer 		wmoNum_use(MAX_RAOBS)
      character*6       staName_use(MAX_RAOBS)
      integer 	timeSyn(MAX_RAOBS),
     1			wmoStaNum(MAX_RAOBS),
     1			timeRel(MAX_RAOBS),
     1			numSigT(MAX_RAOBS),
     1			numSigW(MAX_RAOBS) 
      character*6       staName(MAX_RAOBS) 
      real            staLat(MAX_RAOBS), 
     1			staLon(MAX_RAOBS), 
     1                  staElev(MAX_RAOBS)
      character*9	a9_time(MAX_RAOBS)

C     temp variables
      integer 	numT(MAX_RAOBS),
     1                  timeLapsT(maxT,MAX_RAOBS),
     1			fileAvailTHR
      character*1	typeT(maxT,MAX_RAOBS)
      real            tSigT(maxT,MAX_RAOBS),
     1                  tdSigT(maxT,MAX_RAOBS),
     1                  htSigT(maxT,MAX_RAOBS),
     1                  prSigT(maxT,MAX_RAOBS),
     1			prIT(maxT,MAX_RAOBS), 
     1                  tIT(maxT,MAX_RAOBS),
     1                  tdIT(maxT,MAX_RAOBS),
     1                  prPT(maxT,MAX_RAOBS), 
     1                  tPT(maxT,MAX_RAOBS),
     1                  tdPT(maxT,MAX_RAOBS),
     1                  riT(maxT,MAX_RAOBS),
     1                  rjT(maxT,MAX_RAOBS), 
     1                  rkT(maxT,MAX_RAOBS),
     1                  latT(maxT,MAX_RAOBS), 
     1                  lonT(maxT,MAX_RAOBS), 
     1                  htT(maxT,MAX_RAOBS) 

C     wind variables
      integer 	numW(MAX_RAOBS),
     1                  timeLapsW(maxW,MAX_RAOBS),
     1                  fileAvailUV
      character*1	typeW(maxW,MAX_RAOBS)
      real            htSigW(maxW,MAX_RAOBS), 
     1                  wsSigW(maxW,MAX_RAOBS), 
     1                  wdSigW(maxW,MAX_RAOBS), 
     1                  uIW(maxW,MAX_RAOBS), 
     1                  vIW(maxW,MAX_RAOBS),
     1                  uPW(maxW,MAX_RAOBS), 
     1                  vPW(maxW,MAX_RAOBS),
     1                  riW(maxW,MAX_RAOBS),
     1                  rjW(maxW,MAX_RAOBS), 
     1                  rkW(maxW,MAX_RAOBS),
     1                  latW(maxW,MAX_RAOBS), 
     1                  lonW(maxW,MAX_RAOBS), 
     1                  htW(maxW,MAX_RAOBS) 

C     data written out
      integer         nHts(MAX_RAOBS),
     1                  timeLaps(MAX_HTS,MAX_RAOBS)
      character*1	type(MAX_HTS,MAX_RAOBS)
      real            ri(MAX_HTS,MAX_RAOBS),
     1                  rj(MAX_HTS,MAX_RAOBS), 
     1                  rk(MAX_HTS,MAX_RAOBS),
     1                  lat(MAX_HTS,MAX_RAOBS), 
     1                  lon(MAX_HTS,MAX_RAOBS), 
     1                  hts(MAX_HTS,MAX_RAOBS),
     1			lapsTime(MAX_HTS,MAX_RAOBS),
     1                  uP(MAX_HTS,MAX_RAOBS), 
     1                  uI(MAX_HTS,MAX_RAOBS), 
     1                  vP(MAX_HTS,MAX_RAOBS), 
     1                  vI(MAX_HTS,MAX_RAOBS),
     1                  tP(MAX_HTS,MAX_RAOBS),
     1                  tI(MAX_HTS,MAX_RAOBS), 
     1                  tdP(MAX_HTS,MAX_RAOBS),
     1                  tdI(MAX_HTS,MAX_RAOBS),
     1                  prP(MAX_HTS,MAX_RAOBS), 
     1                  prI(MAX_HTS,MAX_RAOBS)

C     Laps data read in
      real            uLapsGS(ni,nj,nk), vLapsGS(ni,nj,nk),
     1                  tLapsGS(ni,nj,nk), rhLapsGS(ni,nj,nk),
     1                  htLapsGS(ni,nj,nk),uLapsGP(ni,nj,nk),
     1                  vLapsGP(ni,nj,nk),tLapsGP(ni,nj,nk),
     1                  rhLapsGP(ni,nj,nk), htLapsGP(ni,nj,nk),
     1                  htLgaGS(ni,nj,nk), htLgaGP(ni,nj,nk)
       
      integer		writeT,writeW
C
C     BEGIN
C
      istatus = 1   !assume good return
      writeW = 0
      writeT = 0

C     Set up nHts(MAX_RAOBS)
      do i = 1, MAX_RAOBS
        nHts(i) = MAX_HTS
      enddo

C     Read verify_raob.nl

      lun = 10
      call s_len(nl_dir,dir_len)
      filename = nl_dir(1:dir_len)//'/verif_raob.txt'
      open(lun,file=filename,status='old',err=900)

      l_eof = .false.
      n_raobs_use = 0
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
      n_raobs_use = i - 1

      all_raobs = 0	!set to 1 if use all raobs in file

C     See if wmoNum_use(1) .eq. -1...if so, use all raobs in file
      if (wmoNum_use(1) .eq. -1) then
        all_raobs = 1
      endif

C     Read raob file "raob_fname" and return requested raob data
C     with Man/SigT mingled by pressure and Man/SigW mingled by height
      call get_raob_data_a(ni, nj, i4time,i4time_raob_earliest,
     1                     i4time_raob_latest,raob_fname,MAX_RAOBS,
     1                     maxM,maxT, maxW,lats,lons,timeSyn,timeRel, 
     1                     numSigT, numSigW, wmoStaNum,staname,typeW,
     1                     typeT,prSigT, tSigT, tdSigT, htSigT, htSigW,
     1                     wdSigW,wsSigW, staLat, staLon, staElev,
     1                     max_ht_m_proc, min_pres_mb_proc,
     1                     nRaobs, n_raobs_avail,verif_missing_data,
     1                     raob_missing_data, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' Unable to read raob file: ',raob_fname
        istatus = 0
        return
      endif

      write(6,*)
      write(6,*) 'Number of raobs available: ',n_raobs_avail
      write(6,*)

      if (n_raobs_avail .gt. 0) then
   
C       Set which raobs to pull LAPS/model data from
        if (all_raobs .ne. 1) then   !set use_raob(i) to 1 if raob is on list
          do i = 1, n_raobs_use
            do j = 1, nRaobs
              if ((wmoNum_use(i) .eq. wmoStaNum(j)) .or.
     1            (staName_use(i) .eq. staName(j))) then
                use_raob(j) = 1    !set to 1 if current raob is on list to use
              else
                use_raob(j) = 0  
              endif
            enddo
          enddo
        else
          do i = 1, nRaobs
            use_raob(i) = 1    !set to 1 if current raob is on list to use
          enddo
        endif

C       For i4time Read LAPS LW3 file for U and V/ LT1 file for T and HT
C       and LH3 file for rh
C       statusL(1,x) = Syn  (2,x) = Prev
C       statusL(x,1)=u (x,2)=v (x,3)=t (x,4)=ht (x,5)=rh (x,6)=LGAht

        call make_fnam_lp(i4time,a9_time,istatus)

        call get_laps(ni,nj,nk,model_dir,i4time_sys,a9_time,
     1              uLapsGS,vLapsGS,tLapsGS,rhLapsGs,htLapsGS,
     1                uLapsGP,vLapsGP,tLapsGP,rhLapsGP,
     1                htLapsGP,htLgaGS,htLgaGP,balance,
     1                laps_levels_mb,statusL)

        fileAvailUV = 0
        if ((statusL(1,1) .eq. 1) .and. (statusL(1,2) .eq. 1)) then
          fileAvailUV = fileAvailUV + 1
        elseif ((statusL(2,1) .eq. 1) .and.
     1          (statusL(2,2) .eq. 1)) then
          fileAvailUV = fileAvailUV + 2
        else
          istatus = 0
          write(6,*) ' Missing U/V data ',statusL(1,1),
     1    statusL(1,2),' 11Z ',statusL(2,1), statusL(2,2)
        endif

        fileAvailTHR = 0
        if ((statusL(1,3) .eq. 1) .and. (statusL(1,4) .eq. 1) .and.
     1      (statusL(1,5) .eq. 1)) then
          fileAvailTHR = fileAvailTHR + 1
        elseif ((statusL(2,3) .eq. 1) .and.
     1          (statusL(2,4) .eq. 1) .and.
     1          (statusL(2,5) .eq. 1)) then
          fileAvailTHR = fileAvailTHR + 2
        else
          istatus = 0
          write(6,*) ' Missing T/HT/Td data ',statusL(1,1),
     1    statusL(1,2),' 11Z ',statusL(2,1), statusL(2,2)
        endif

        if (fileAvailUV .gt. 0) then !need UV for both..if not there, can't run

C         Loop through use_raob and pull model data if use_raob(i) .eq. 1
          do i = 1, nRaobs
            if (use_raob(i) .eq. 1) then

C             Ascend balloon for W and return ri,rj,rk,latW,lonW,htW, uIW,vIW,
C                                           timeLapsW,istatus
              if (fileAvailUV .gt. 0) then

                call ascend_w(timeSyn(i), timeRel(i), numSigW(i),
     1                      numW(i),fileAvailUV,wmoStaNum(i),
     1                      staLat(i), staLon(i), staElev(i), maxW,
     1                      max_ht_m_proc, typeW,
     1                      MAX_RAOBS, ni, nj, nk, i, statusL,
     1                      htSigW, wdSigW, wsSigW,  !remember (maxW,MAX_RAOBS)
     1                      uLapsGS,vLapsGS,uLapsGP,vLapsGP,  !(ni,nj,nk)
     1                      htLapsGS, htLapsGP, raob_missing_data,
     1                      verif_missing_data,
     1                      lats,lons,laps_levels_pa,    !(ni,nj,nk)
!................................variables below returned...................
     1                      riW,rjW,rkW,latW,lonW,htW, uIW,vIW,
     1                      uPW,vPW,timeLapsW,istatus)

                if (istatus .ne. 1) then
                  write(6,*) 'Wind data for raob ',wmoStaNum(i),
     1                     ' not available.'
           
                  numW(i) = 0
                endif
              else
                numW(i) = 0
              endif

              if (fileAvailTHR .gt. 0) then

                call ascend_t(timeSyn(i), timeRel(i), numSigT(i),
     1                        numT(i),fileAvailUV,fileAvailTHR, 
     1                        wmoStaNum(i),staLat(i), staLon(i), 
     1                        staElev(i), typeT,
     1                        maxT, maxW,max_ht_m_proc,
     1                        MAX_RAOBS, ni, nj, nk, i, statusL,
     1                        prSigT, tSigT, tdSigT,htSigT,  !remember (maxT,MAX_RAOBS)
     1                        tLapsGS,rhLapsGS,htLapsGS,tLapsGP,  !(ni,nj,nk)
     1                        rhLapsGP, htLapsGP, 
     1                        htLgaGS, htLgaGP,raob_missing_data,
     1                        verif_missing_data, lats,lons, !(ni,nj,nk)
     1                        laps_levels_pa, numSigW(i),
     1                        htSigW, wdSigW, wsSigW, !remember (maxW,MAX_RAOBS)
!................................variables below returned...................
     1                        riT,rjT,rkT,latT,lonT,htT,prIT,tIT,
     1                        tdIt,prPT,tPT,tdPT,timeLapsT,istatus)

                if (istatus .ne. 1) then
                  numT(i) = 0
                  write(6,*) 'Temp data for raob ',wmoStaNum(i),
     1                       ' not available.'
                endif
              else
                numT(i) = 0
              endif

              if ((numT(i) .eq. 0) .and. (numW(i) .eq. 0)) then
                use_raob(i) = 0
              else
                if (numT(i) .ne. 0) writeT = writeT + 1
                if (numW(i) .ne. 0) writeW = writeW + 1
              endif

            endif
          enddo

C         re-calc n_raobs_use based on interpolation errors
          n_raobs_use = 0
          do i = 1, nRaobs
            if (use_raob(i) .eq. 1) n_raobs_use = n_raobs_use + 1
          enddo

C         merge T and W data
c         call mergeTW(MAX_HTS,MAX_RAOBS,maxW,maxT,numW,numT,
c    1               nRaobs,use_raob,n_raobs_use,
c    1               ri,rj,rk,lat,lon,hts,type,uP,uI,vP,vI,
c    1               tP,tI,tdP,tdI,prP,prI,timeLaps,nHts,
c    1               riT,rjT,rkT,latT,lonT,htT,typeT,
c    1               riW,rjW,rkW,latW,lonW,htW,typeW, 
c    1               uIW,vIW,uPW,vPW,timeLapsW,
c    1               prIT,tIT,tdIT,prPT,tPT,tdPT,timeLapsT,
c    1               verif_missing_data, raob_missing_data,
c    1               istatus)


C         write output files

          if (n_raobs_use .gt. 0) then
            call write_verif_raob(output_fname, MAX_RAOBS,maxW,
     1                          maxT,numW,numT,nRaobs,
     1                          writeT, writeW,
     1                          n_raobs_use,use_raob, wmoStaNum,
     1                          staName,staLat,staLon, staElev,
     1                          a9_time,timeSyn,timeRel,
     1                          riT,rjT,rkT,latT,lonT,htT,
     1                          riW,rjW,rkW,latW,lonW,htW, 
     1                          uIW,vIW,uPW,vPW,timeLapsW,
     1                          prIT,tIT,tdIT,prPT,tPT,tdPT,
     1                          timeLapsT, istatus)

            if (istatus .ne. 1) then
              write(6,*) ' Error writing raob verif file: ',
     1                    output_fname
            endif
          else
            write(6,*) 'No raobs to output'
          endif

        else !no model data files to process
          write(6,*) 'no UV data available...cannot verify raobs'
        endif
        goto 999

      else
        write(6,*)'No raobs to output' 
        goto 999
      endif

900   print*,'error opening namelist file ', filename
      goto 999

901   print*,'error reading namelist file ', filename
      goto 999


999   continue

      return
      end
!1.........................................................................
      subroutine write_verif_raob(output_fname, MAX_RAOBS, maxW,
     1                          maxT,numW,numT,nRaobs,
     1                          writeT,writeW,
     1                          n_raobs_use,use_raob, wmoStaNum,
     1                          staName,staLat,staLon,staElev,
     1                          a9_time,timeSyn,timeRel,
     1                          riT,rjT,rkT,latT,lonT,htT,
     1                          riW,rjW,rkW,latW,lonW,htW, 
     1                          uIW,vIW,uPW,vPW,timeLapsW,
     1                          prIT,tIT,tdIT,prPT,tPT,tdPT,
     1                          timeLapsT,istatus)

      implicit none

      character*256     output_fname	!path and name of output file
      integer		MAX_RAOBS, maxW,maxT
      integer           nRaobs, n_raobs_use, use_raob(MAX_RAOBS)
      integer		writeT,writeW,wmoStaNum(MAX_RAOBS) 
      character*6       staName(MAX_RAOBS) 
      real            staLat(MAX_RAOBS), staLon(MAX_RAOBS), 
     1                  staElev(MAX_RAOBS)
      character*9	a9_time(MAX_RAOBS)
      integer		timeSyn(MAX_RAOBS),
     1         		timeRel(MAX_RAOBS)

C     temp variables
      integer         numT(MAX_RAOBS),
     1			timeLapsT(maxT,MAX_RAOBS)
      real            riT(maxT,MAX_RAOBS),
     1                  rjT(maxT,MAX_RAOBS),
     1                  rkT(maxT,MAX_RAOBS),
     1                  latT(maxT,MAX_RAOBS),
     1                  lonT(maxT,MAX_RAOBS),
     1                  htT(maxT,MAX_RAOBS),
     1                  prIT(maxT,MAX_RAOBS),
     1                  tIT(maxT,MAX_RAOBS),
     1                  tdIT(maxT,MAX_RAOBS),
     1                  prPT(maxT,MAX_RAOBS),
     1                  tPT(maxT,MAX_RAOBS),
     1                  tdPT(maxT,MAX_RAOBS)

C     wind variables
      integer         numW(MAX_RAOBS),
     1 			timeLapsW(maxW,MAX_RAOBS)
      real            riW(maxW,MAX_RAOBS),
     1                  rjW(maxW,MAX_RAOBS),
     1                  rkW(maxW,MAX_RAOBS),
     1                  latW(maxW,MAX_RAOBS),
     1                  lonW(maxW,MAX_RAOBS),
     1                  htW(maxW,MAX_RAOBS),
     1                  uIW(maxW,MAX_RAOBS),
     1                  vIW(maxW,MAX_RAOBS),
     1                  uPW(maxW,MAX_RAOBS),
     1                  vPW(maxW,MAX_RAOBS)

      integer		istatus

C     local variables
      integer		i, j
      character*256     output_fnameW   !path/name of wind output file
      character*256     output_fnameT   !path/name of temp output file
      integer		output_len

C
C     BEGIN
C
      istatus = 1   !assume good return

      call s_len(output_fname,output_len)
      output_fnameT = output_fname(1:output_len)//'T'

      if (writeT .ne. 0) then
        write(6,*) output_fnameT

C       open output_fname
        open(1,file=output_fnameT,status='unknown',err=98)
        go to 99

98      write(6,*)' Error opening temp verif file: ',output_fnameT
        istatus = 0
        return

99      continue

C 100   writes: n_raobs_use
100     format(i3)
C 101   writes: nHts,a9_time,wmoStaNum,staName,staLat,staLon,staElev,timeSyn,timeRel
101     format(i3,1x,a9,1x,i7,1x,a6,1x,f7.3,1x,f8.3,1x,f7.0,2(1x,i12))
C 102   writes: htT,latT,lonT,riT,rjT,rkT,tPT,tIT,tdPT,tdIT,prPT,prIT,timeLapsT
102     format(f7.1,1x,f6.2,1x,f8.2,1x,3(f8.3,1x),6(f6.1,1x),i12)


C       write number of raob into file
        write(1,100) writeT
      
C       write temp data out for each raob where use_raob(i) .eq. 1
        do i = 1, nRaobs
          if ((use_raob(i) .eq. 1).and.(numT(i).gt.0)) then
C           write header info
            write(1,101)numT(i),a9_time(i),wmoStaNum(i),
     1                 staName(i)(1:6),
     1                 staLat(i),staLon(i),staElev(i),
     1                 timeSyn(i),timeRel(i)


            write(1,104)
            do j = 1, numT(i)
              write(1,102) htT(j,i),latT(j,i),lonT(j,i),
     1                   riT(j,i),rjT(j,i),rkT(j,i),tPT(j,i),
     1                   tIT(j,i),tdPT(j,i),tdIT(j,i),
     1                   prPT(j,i),prIT(j,i),timeLapsT(j,i)
            enddo
          endif
        enddo

        close(1)
      else
        write(6,*) 'No temp data written: ',output_fnameT
      endif
104     format(2x,'hgt',5x,'lat',5x,'lon',6x,'ri',7x,'rj',7x,'rk',
     +7x,'tPT',4x,'tIT',3x,'tdPT',3x,'tdIT',3x,'prPT',3x,'prIT'
     +,3x,'i4timelapsT')

      output_fnameW = output_fname(1:output_len)//'W'
      if (writeW .gt. 0) then
        write(6,*) output_fnameW

C       open output_fname
        open(1,file=output_fnameW,status='unknown',err=198)
        go to 199

198     write(6,*)' Error opening wind verif file: ',output_fnameW
        istatus = 0
        return

199     continue

C 103   writes: htW,latW,lonW,riW,rjW,rkW,uPW,uIW,vPW,vIW,timeLapsW
103     format(f7.1,1x,f6.2,1x,f8.2,1x,3(f8.3,1x),4(f6.1,1x),i12)

C       write number of raob into file
        write(1,100) n_raobs_use
      
C       write wind data out for each raob where use_raob(i) .eq. 1
        do i = 1, nRaobs
          if ((use_raob(i) .eq. 1).and.(numW(i).gt.0)) then
C           write header info
            write(1,101)numW(i),a9_time(i),wmoStaNum(i),
     1                 staName(i)(1:6),
     1                 staLat(i),staLon(i),staElev(i),
     1                 timeSyn(i),timeRel(i)

            write(1,105)
            do j = 1, numW(i)
            write(1,103) htW(j,i),latW(j,i),lonW(j,i),riW(j,i),
     1                   rjW(j,i),rkW(j,i),uPW(j,i),
     1                   uIW(j,i),vPW(j,i),vIW(j,i),timeLapsW(j,i)
            enddo
          endif
        enddo

        close(1)
      else
        write(6,*) 'No wind data to write out: ',output_fnameW
      endif

105   format(2x,'hgt',5x,'lat',5x,'lon',6x,'ri',7x,'rj',7x,'rk',
     +7x,'uPW',4x,'uIW',4x,'vPW',4x,'vIW',3x,'i4timelapsW')

      return
      end
!2.........................................................................
