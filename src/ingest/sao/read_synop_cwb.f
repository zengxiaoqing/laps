      subroutine read_synop_cwb ( filename, maxSkyCvr, maxobs,
     ~                          i4time_sys, path_to_local,            
     ~                          altm, stnTp, td, tdTths, elev,
     ~                          lats, lons, t24max, t24min,
     ~                          pcp1hr, pcp24hr, pcp3hr, pcp6hr, 
     ~                          prsWth, p, pc, pcc, rptTp, rh,
     ~                          mslp, skyCvr, skyLyrBs, snowCvr, sr, st,       
     ~                          stname, tTths, t, timeObs, vis,
     ~                          dd, wgdd, ff, wgff, wmoId, badflag,
     ~                          num, istatusSynop )

      integer, parameter :: maxSynop = 150
      integer, parameter :: maxMso =    50

      character*(*)  filename, path_to_local
      character(25)  prsWth(maxobs)
      character(13)  cvt_i4time_wfo_fname13,a13time_eat
      character(9)   a9time
      character(8)   skyCvr(maxSkyCvr,maxobs)
      character(6)   rptTp(maxobs), stnTp(maxobs)
      character(5)   stname(maxobs), stnNo(maxobs)

      integer  pcc(maxobs), wmoId(maxobs), flag

      real  altm(maxobs), td(maxobs), tdTths(maxobs)
      real  elev(maxobs), lats(maxobs), lons(maxobs)
      real  t24max(maxobs), t24min(maxobs)
      real  rh(maxobs), pcp1hr(maxobs), pcp24hr(maxobs)
      real  pcp3hr(maxobs), pcp6hr(maxobs)
      real  p(maxobs), pc(maxobs), mslp(maxobs)
      real  skyLyrBs(maxSkyCvr,maxobs), snowCvr(maxobs)
      real  tTths(maxobs), t(maxobs)
      real  vis(maxobs), dd(maxobs), wgdd(maxobs), wgff(maxobs)
      real  ff(maxobs), sr(maxobs), st(maxobs)

      double precision  timeObs(maxobs), timeObsMso(maxMso)

      character(6)   rptTpMso(maxobs), stnTpMso(maxobs)
      character(5)   stnNoMso(maxMso)
      integer  pccMso(maxMso)

      real  latsMso(maxMso), lonsMso(maxMso)
      real  elevMso(maxMso), tMso(maxMso), t24maxMso(maxMso)
      real  t24minMso(maxMso), tdMso(maxMso), rhMso(maxMso)
      real  pcp1hrMso(maxMso), pcp3hrMso(maxMso), pcp6hrMso(maxMso)
      real  pcp24hrMso(maxMso), ddMso(maxMso), ffMso(maxMso)
      real  wgddMso(maxMso), wgffMso(maxMso), pMso(maxMso)
      real  mslpMso(maxMso), pcMso(maxMso), srMso(maxMso), stMso(maxMso)       

      rptTp   = 'SYNOP'
      stnTp   = 'UNK'
      skyCvr  = '        '
      stname  = 'UNK'
      stnNo   = '     '
      wmoId   = ibadflag
      pcc     = ibadflag
      timeObs = badflag
      lats    = badflag
      lons    = badflag
      elev    = badflag
      t       = badflag
      tTths   = badflag
      t24max  = badflag
      t24min  = badflag
      td      = badflag
      tdTths  = badflag
      rh      = badflag
      pcp1hr  = badflag
      pcp3hr  = badflag
      pcp6hr  = badflag
      pcp24hr = badflag
      ff      = badflag
      dd      = badflag
      wgff    = badflag
      wgdd    = badflag
      p       = badflag
      mslp    = badflag
      pc      = badflag
      altm    = badflag
      skyLyrBs= badflag
      vis     = badflag
      snowCvr = badflag
      sr      = badflag
      st      = badflag
      
      istatusSynop= 0
      istatusMso=   0
      num=          0
      numSynop=     0
      numMso=       0

      nq= maxSynop

      call read_synop_cwb_sub ( filename, maxSkyCvr, maxSynop,
     ~     altm(1:nq), stnTp(1:nq), td(1:nq), tdTths(1:nq), elev(1:nq),
     ~     lats(1:nq), lons(1:nq), t24max(1:nq), t24min(1:nq),
     ~     pcp1hr(1:nq), pcp24hr(1:nq), pcp3hr(1:nq), pcp6hr(1:nq), 
     ~     prsWth(1:nq), pc(1:nq), pcc(1:nq), rptTp(1:nq), mslp(1:nq),
     ~     skyCvr(1:maxSkyCvr,1:nq), skyLyrBs(1:maxSkyCvr,1:nq),
     ~     snowCvr(1:nq), stname(1:nq), tTths(1:nq), t(1:nq),
     ~     timeObs(1:nq), vis(1:nq), dd(1:nq), wgff(1:nq), ff(1:nq),
     ~     wmoId(1:nq), badflag, numSynop, istatusSynop )

      do i= 1,numSynop
         write(stnNo(i),'(i5)') wmoId(i)      
c        write(*,*) i,stnNo(i),rptTp(i),stnTp(i), timeObs(i)
         write(*,*) i,stnNo(i),stname(i)
      enddo

      np= numSynop +1
      nq= numSynop +maxMso

      call read_meso_cwb ( path_to_local, maxMso, badflag, ibadflag, 
     ~                     i4time_sys, timeObsMso, rptTpMso, stnTpMso,
     ~                     stnNoMso, latsMso, lonsMso, elevMso,
     ~                     tMso, t24maxMso, t24minMso, tdMso, rhMso, 
     ~                     pcp1hrMso, pcp3hrMso, pcp6hrMso, pcp24hrMso,
     ~                     ddMso, ffMso, wgddMso, wgffMso, pMso, 
     ~                     mslpMso, pccMso, pcMso, srMso, stMso, 
     ~                     numMso, istatusMso )

c                    combine synop data and mesonet data 
      k= numSynop
      do i= 1,numMso
         flag= 0

c write(*,*) i,stnNoMso(i),rptTpMso(i),stnTpMso(i),timeObsMso(i)
         do j= 1,numSynop
            if ( stnNoMso(i) .eq. stnNo(j) ) then
               rh(j)=   rhMso(i)
               dd(j)=   ddMso(i)
               ff(j)=   ffMso(i)
               wgdd(j)= wgddMso(i)
               wgff(j)= wgffMso(i)
               p(j)=    pMso(i)
               sr(j)=   srMso(i)
               st(j)=   stMso(i)

               flag= 1
               exit
            endif
         enddo 

         if ( flag /= 1 ) then
            k= k +1

            timeObs(k)= timeObsMso(i)
            rptTp(k)=   rptTpMso(i)
            stnTp(k)=   stnTpMso(i)
            stnNo(k)=   stnNoMso(i)
            wmoId(k)=   ibadflag
            td(k)=      tdMso(i)
            tdTths(k)=  tdMso(i)
            elev(k)=    elevMso(i)
            lats(k)=    latsMso(i)
            lons(k)=    lonsMso(i)
            t24max(k)=  t24maxMso(i)
            t24min(k)=  t24minMso(i)
            rh(k)=      rhMso(i)
            rh(k)=      rhMso(i)
            pcp1hr(k)=  pcp1hrMso(i)
            pcp24hr(k)= pcp24hrMso(i)
            pcp3hr(k)=  pcp3hrMso(i)
            pcp6hr(k)=  pcp6hrMso(i)
            p(k)=       pMso(i)
            pc(k)=      pcMso(i)
            pcc(k)=     pccMso(i)
            mslp(k)=    mslpMso(i)
            t(k)=       tMso(i)
            dd(k)=      ddMso(i)
            wgdd(k)=    wgddMso(i)
            ff(k)=      ffMso(i)
            wgff(k)=    wgffMso(i)
            sr(k)=      srMso(i)
            st(k)=      stMso(i)
         endif
      enddo

      do i= 1,k
      write(*,*) i,stnNo(i),timeObs(i),rptTp(i),stnTp(i),
     ~           stname(i),wmoId(i), stnNo(i)
      enddo
      num = k

      end



      subroutine read_synop_cwb_sub (
     ~     filename, maxSkyCover, recNum, altimeter,
     ~     autoStationType, dewpoint, dpFromTenths, elevation,
     ~     latitude, longitude, maxTemp24Hour, minTemp24Hour,
     ~     precip1Hour, precip24Hour, precip3Hour, precip6Hour,
     ~     presWeather, pressChange3Hour, pressChangeChar,
     ~     reportType, seaLevelPress, skyCover, skyLayerBase,
     ~     snowCover, stationName, tempFromTenths, temperature,
     ~     timeObs, visibility, windDir, windGust, windSpeed, wmoId,
     ~     badflag, staNum, istatus )

      integer  maxSkyCover, recNum

      character*(*)  filename
      character(6)   autoStationType(recNum)
      character(25)  presWeather(recNum)
      character(6)   reportType(recNum)
      character(8)   skyCover(maxSkyCover,recNum)
      character(5)   stationName(recNum)

      integer  pressChangeChar(recNum), wmoId(recNum)

      real  altimeter(recNum), dewpoint(recNum), dpFromTenths(recNum)
      real  elevation(recNum), latitude(recNum), longitude(recNum)
      real  maxTemp24Hour(recNum), minTemp24Hour(recNum)
      real  precip1Hour(recNum), precip24Hour(recNum)
      real  precip3Hour(recNum), precip6Hour(recNum)
      real  pressChange3Hour(recNum), seaLevelPress(recNum)
      real  skyLayerBase(maxSkyCover,recNum), snowCover(recNum)
      real  tempFromTenths(recNum), temperature(recNum)
      real  visibility(recNum), windDir(recNum), windGust(recNum)
      real  windSpeed(recNum)
      real  badflag

      double precision  timeObs(recNum)

      character(3)   reportFlag(recNum)
      character(2)   yy(recNum), mo(recNum), dd(recNum)
      character(2)   hh(recNum), mn(recNum)
      character(10)  time(recNum)
      character(9)   a10_to_a9

      integer  windQua(recNum), seaLevelPressQua(recNum)
      integer  temperatureQua(recNum), dewpointQua(recNum)
      integer  pressChange3HourQua(recNum)
      integer  dupliStation(9), staNum, dummy

      real  tempDewDiff(recNum), lowestCloudHeight(0:10)

      istatus= 0
      staNum= 0

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      data  dupliStation / 58968, 58974, 59158, 59358, 59559,
     ~           	   59562, 59567, 59792, 59997 /
      data  lowestCloudHeight / 0., 50., 100., 200., 300., 600., 1000.,
     ~                          1500., 2000., 2500., 0. /

c      ------   give initial values to avoid data stack problem  ------
      do 5 j= 1,recNum
      do 5 i= 1,maxSkyCover
         skyCover(i,j)= "   "
5        skyLayerBase(i,j)= badflag

      do j= 1,recNum
         read (1,20,end=99,err=9) reportFlag(j), wmoId(j),
     ~                   elevation(j), latitude(j), longitude(j),
     ~                   yy(j), mo(j), dd(j), hh(j), mn(j)
         read (1,30,end=9,err=9) windDir(j), windSpeed(j), windQua(j),
     ~                   visibility(j), presWeather(j),
     ~                   seaLevelPress(j), seaLevelPressQua(j),
     ~                   temperature(j), temperatureQua(j),
     ~                   skyCover(1,j), skyLayerBase(1,j)
         read (1,40,end=9,err=9) tempDewDiff(j), dewpointQua(j),       
     ~                   pressChangeChar(j), pressChange3Hour(j),
     ~                   pressChange3HourQua(j), precip3Hour(j),
     ~                   maxTemp24Hour(j), minTemp24Hour(j), windGust(j) 
         read (1,50,end=9,err=9) skyCover(2,j), skyLayerBase(2,j),
     ~                   skyCover(3,j), skyLayerBase(3,j), 
     ~                   precip24Hour(j)
         read (1,*)

         if ( reportFlag(j) .ne. '*31' )  then
            write (6,*) 'read synop data heading error'
            go to 1000
         endif
	 go to 10

9        write (6,*) ' Reading error of synop data'
         write (6,*) reportFlag(j), wmoId(j), elevation(j),
     ~               latitude(j), longitude(j), ' ',
     ~               yy(j), mo(j), dd(j), hh(j), mn(j), ' ', timeobs(j)
         write (6,*) windDir(j), windSpeed(j), windQua(j),
     ~               visibility(j), presWeather(j),
     ~               seaLevelPress(j), seaLevelPressQua(j),
     ~               temperature(j), temperatureQua(j),
     ~               skyCover(1,j), skyLayerBase(1,j)
         write (6,*) dewpoint(j), dewpointQua(j),
     ~               pressChangeChar(j), pressChange3Hour(j),
     ~               pressChange3HourQua(j), precip3Hour(j),
     ~               maxTemp24Hour(j), minTemp24Hour(j), windGust(j) 
         write (6,*) skyCover(2,j), skyLayerBase(2,j), 
     ~               skyCover(3,j), skyLayerBase(3,j), precip24Hour(j)

10       staNum= staNum +1
      enddo

20    format ( a3, i5, f4.0, 2f5.2, 2x, 5a2 )
30    format ( 2x, 2f3.0, i1, f2.0, a2, 3x, f5.1, i1, f4.1, i1, a2, 2x,
     ~         f2.0 )
40    format ( f3.1, i1, 1x, i2, f3.1, i1, 3(1x, f4.1), 8x, f3.0 )
50    format ( 2(a2, 2x, f2.0), 8x, f4.1 )

c     --- eliminate duplicate data coming from international broadcast ---
99    do 100 k= 1,9
      do 100 j= 1,staNum
         if ( wmoId(j) .eq. dupliStation(k) )  then
            staNum= staNum -1
            do i= j,staNum
               reportFlag(i)= reportFlag(i+1)
	       wmoId(i)= wmoId(i+1)
	       elevation(i)= elevation(i+1)
	       latitude(i)= latitude(i+1)
	       longitude(i)= longitude(i+1)
	       yy(i)= yy(i+1)
	       mo(i)= mo(i+1) 
	       dd(i)= dd(i+1)
	       hh(i)= hh(i+1)
	       mn(i)= mn(i+1)
	       windDir(i)= windDir(i+1)
	       windSpeed(i)= windSpeed(i+1)
	       windQua(i)= windQua(i+1)
	       visibility(i)= visibility(i+1)
	       presWeather(i)= presWeather(i+1)
	       seaLevelPress(i)= seaLevelPress(i+1)
	       seaLevelPressQua(i)= seaLevelPressQua(i+1)
	       temperature(i)= temperature(i+1)
	       temperatureQua(i)= temperatureQua(i+1)
	       skyCover(1,i)= skyCover(1,i+1)
	       skyLayerBase(1,i)= skyLayerBase(1,i+1)
	       tempDewDiff(i)= tempDewDiff(i+1)
	       dewpointQua(i)= dewpointQua(i+1)
	       pressChangeChar(i)= pressChangeChar(i+1)
	       pressChange3Hour(i)= pressChange3Hour(i+1)
	       pressChange3HourQua(i)= pressChange3HourQua(i+1)
	       precip3Hour(i)= precip3Hour(i+1)
	       maxTemp24Hour(i)= maxTemp24Hour(i+1)
	       minTemp24Hour(i)= minTemp24Hour(i+1)
	       windGust(i)= windGust(i+1)
               skyCover(2,i)= skyCover(2,i+1)
	       skyLayerBase(2,i)= skyLayerBase(2,i+1)
	       skyCover(3,i)= skyCover(3,i+1)
	       skyLayerBase(3,i)= skyLayerBase(3,i+1)
	       precip24Hour(i)= precip24Hour(i+1)
            enddo
         endif
100   continue

c      ----------       examine data quality and change units       ---------
      do j= 1,staNum
         if ( windQua(j) .ne. 1 )  then
            windDir(j)= badflag
            windSpeed(j)= badflag
         endif

         if ( windGust(j) .eq. -99. )  windGust(j)= badflag
         if ( elevation(j) .eq. -999. )  elevation(j)= badflag
         if ( pressChangeChar(j).eq.-9 ) pressChangeChar(j)=int(badflag)
         if ( presWeather(j) .eq. '-9' )  presWeather(j)= 'UNK'

         if ( seaLevelPressQua(j) .eq. 1 )  then
               seaLevelPress(j)= seaLevelPress(j) *100.   ! millibar -> pascal
            else
               seaLevelPress(j)= badflag
         endif

         if ( pressChange3HourQua(j) .eq. 1 )  then
            pressChange3Hour(j)= pressChange3Hour(j) *100. ! millibar -> pascal
          else
            pressChange3Hour(j)= badflag
            pressChangeChar(j)= int(badflag)
         endif

         if ( temperatureQua(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15           ! degC -> degK
            else
               temperature(j)= badflag
         endif

         tempFromTenths(j)= temperature(j)

         if ( dewpointQua(j) .eq. 1 )  then
               dewpoint(j)= temperature(j) -tempDewDiff(j)        ! unit: degK
            else
               dewpoint(j)= badflag
         endif

         dpFromTenths(j)= dewpoint(j)

         if ( maxTemp24Hour(j) .eq. -99.9 )  then
               maxTemp24Hour(j)= badflag
            else
               maxTemp24Hour(j)= maxTemp24Hour(j) +273.15       ! degC -> degK
         endif

         if ( minTemp24Hour(j) .eq. -99.9 )  then
               minTemp24Hour(j)= badflag
            else
               minTemp24Hour(j)= minTemp24Hour(j) +273.15       ! degC -> degK
         endif

         if ( precip24Hour(j) .eq. -99.9 )  then
               precip24Hour(j)= badflag
            else
               precip24Hour(j)= precip24Hour(j) *0.001   ! millimeter -> meter
         endif

         if ( precip3Hour(j) .eq. -99.9 )  then
               precip3Hour(j)= badflag
            else
               precip3Hour(j)= precip3Hour(j) *0.001     ! millimeter -> meter
         endif

         if ( yy(j)(1:1) .eq. ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( mo(j)(1:1) .eq. ' ' )  mo(j)= '0'//mo(j)(2:2)
         if ( dd(j)(1:1) .eq. ' ' )  dd(j)= '0'//dd(j)(2:2)
         if ( hh(j)(1:1) .eq. ' ' )  hh(j)= '0'//hh(j)(2:2)
         if ( mn(j)(1:1) .eq. ' ' )  mn(j)= '0'//mn(j)(2:2)
         time(j)= yy(j)//mo(j)//dd(j)//hh(j)//mn(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeObs(j)= dble( i4time )                       ! seconds since 1960
      enddo

c    -------    transform code figure into visibility ( unit: m )  -------
      do j= 1,staNum
         if     ( visibility(j) .eq.  0. )  then
               visibility(j)= 50.
         elseif ( visibility(j) .gt.  0.  .and. 
     ~            visibility(j) .lt. 51. )  then
               visibility(j)= visibility(j) *100.
         elseif ( visibility(j) .gt. 55.  .and. 
     ~            visibility(j) .lt. 81. )  then
               visibility(j)= ( visibility(j) -50. ) *1000.
         elseif ( visibility(j) .gt. 80.  .and. 
     ~            visibility(j) .lt. 90. )  then
               visibility(j)= ( visibility(j) -74. ) *5000.
         elseif ( visibility(j) .eq. 90. )  then
               visibility(j)= 25.
         elseif ( visibility(j) .eq. 91. )  then
               visibility(j)= 50.
         elseif ( visibility(j) .eq. 92. )  then
               visibility(j)= 200.
         elseif ( visibility(j) .eq. 93. )  then
               visibility(j)= 500.
         elseif ( visibility(j) .eq. 94. )  then
               visibility(j)= 1000.
         elseif ( visibility(j) .eq. 95. )  then
               visibility(j)= 2000.
         elseif ( visibility(j) .eq. 96. )  then
               visibility(j)= 4000.
         elseif ( visibility(j) .eq. 97. )  then
               visibility(j)= 10000.
         elseif ( visibility(j) .eq. 98. )  then
               visibility(j)= 20000.
         elseif ( visibility(j) .eq. 99. )  then
               visibility(j)= 50000.
         else
               visibility(j)= badflag
         endif
      enddo

c ---- transform code figure into base of cloud layer indicated (unit: m) ----
      do 200 i= 2,maxSkyCover
      do 200 j= 1,staNum
         if     ( skyLayerBase(i,j) .eq.  0. )  then
               skyLayerBase(i,j)= 15.
         elseif ( skyLayerBase(i,j) .gt.  0.  .and.  
     ~            skyLayerBase(i,j) .lt. 51. )  then
               skyLayerBase(i,j)= skyLayerBase(i,j) *30.
         elseif ( skyLayerBase(i,j) .gt. 55.  .and.  
     ~            skyLayerBase(i,j) .lt. 81. )  then
               skyLayerBase(i,j)= ( skyLayerBase(i,j) -50. ) *300.
         elseif ( skyLayerBase(i,j) .gt. 80.  .and. 
     ~            skyLayerBase(i,j) .lt. 90. )  then
               skyLayerBase(i,j)= ( skyLayerBase(i,j) -74. ) *1500.
         elseif ( skyLayerBase(i,j) .eq. 90. )  then
               skyLayerBase(i,j)= 25.
         elseif ( skyLayerBase(i,j) .eq. 91. )  then
               skyLayerBase(i,j)= 75.
         elseif ( skyLayerBase(i,j) .eq. 92. )  then
               skyLayerBase(i,j)= 150.
         elseif ( skyLayerBase(i,j) .eq. 93. )  then
               skyLayerBase(i,j)= 250.
         elseif ( skyLayerBase(i,j) .eq. 94. )  then
               skyLayerBase(i,j)= 450.
         elseif ( skyLayerBase(i,j) .eq. 95. )  then
               skyLayerBase(i,j)= 800.
         elseif ( skyLayerBase(i,j) .eq. 96. )  then
               skyLayerBase(i,j)= 1250.
         elseif ( skyLayerBase(i,j) .eq. 97. )  then
               skyLayerBase(i,j)= 1750.
         elseif ( skyLayerBase(i,j) .eq. 98. )  then
               skyLayerBase(i,j)= 2250.
         elseif ( skyLayerBase(i,j) .eq. 99. )  then
               skyLayerBase(i,j)= 3000.
         else
               skyLayerBase(i,j)= badflag
         endif
200   continue

      do 300 i= 2,maxSkyCover
      do 300 j= 1,staNum

!        Modified by Steve Albers  5/9/2001
!        if ( skyCover(i,j).eq.'-9' .or. skyLayerBase(i,j).eq.'-9' )
         if ( skyCover(i,j).eq.'-9' .or. skyLayerBase(i,j).eq.-9 )
     ~      go to 300

c          ----- eliminate duplicate skyCovers and skyLayerBases -----
         dummy= int( skyLayerBase(1,j) )
         if ( dummy .gt. 0 )  then
            if ( skyCover(i,j) .eq. skyCover(1,j)  .and.
     ~           skyLayerBase(i,j) .ge. lowestCloudHeight(dummy) )  then
               if ( dummy .ne. 9  .and.
     ~              skyLayerBase(i,j) .gt. lowestCloudHeight(dummy+1) )
     ~            go to 250
               skyCover(1,j)= '   '
	       skyLayerBase(1,j)= badflag
            endif
         endif

c        ----- eliminate unreasonable skyCovers and skyLayerBases -----
250      if ( skyCover(i,j).eq.' 8' .and. skyCover(1,j).eq.' 8' )  then
	    if ( skyLayerBase(i,j) .gt. skyLayerBase(1,j) )  then
               skyCover(i,j)= '   '
               skyLayerBase(i,j)= badflag
            else
	       skyCover(1,j)= '   '
	       skyLayerBase(1,j)= badflag
            endif
         endif
300   continue

c  -----  transform code figure into base of lowest cloud ( unit: m )  -----
      do j= 1,staNum
	 if     ( skyLayerBase(1,j) .eq. 0. )  then
       	       skyLayerBase(1,j)= 25.
         elseif ( skyLayerBase(1,j) .eq. 1. )  then
               skyLayerBase(1,j)= 75.
	 elseif ( skyLayerBase(1,j) .eq. 2. )  then
	       skyLayerBase(1,j)= 150.
         elseif ( skyLayerBase(1,j) .eq. 3. )  then
	       skyLayerBase(1,j)= 250.
	 elseif ( skyLayerBase(1,j) .eq. 4. )  then
	       skyLayerBase(1,j)= 450.
	 elseif ( skyLayerBase(1,j) .eq. 5. )  then
	       skyLayerBase(1,j)= 800.
	 elseif ( skyLayerBase(1,j) .eq. 6. )  then
	       skyLayerBase(1,j)= 1250.
	 elseif ( skyLayerBase(1,j) .eq. 7. )  then
	       skyLayerBase(1,j)= 1750.
	 elseif ( skyLayerBase(1,j) .eq. 8. )  then
	       skyLayerBase(1,j)= 2250.
	 elseif ( skyLayerBase(1,j).eq.9. .and. skyCover(1,j).ne.' 0'
     ~             .and. skyCover(1,j).ne.'-9' )  then
	       skyLayerBase(1,j)= 3000.
	 else
	       skyLayerBase(1,j)= badflag
  	 endif
      enddo
 
c        --- transform code figure of cloud cover into metar format ---
      do 400 i= 1,maxSkyCover
      do 400 j= 1,staNum
         if     ( skyCover(i,j) .eq. ' 0' )  then
            skyCover(i,j)= 'CLR'
         elseif ( skyCover(i,j) .eq. ' 1'  .or.
     ~            skyCover(i,j) .eq. ' 2' )  then
            skyCover(i,j)= 'FEW'
         elseif ( skyCover(i,j) .eq. ' 3'  .or.
     ~            skyCover(i,j) .eq. ' 4' )  then
            skyCover(i,j)= 'SCT'
         elseif ( skyCover(i,j) .eq. ' 5'  .or.
     ~            skyCover(i,j) .eq. ' 6'  .or.
     ~            skyCover(i,j) .eq. ' 7' )  then
            skyCover(i,j)= 'BKN'
         elseif ( skyCover(i,j) .eq. ' 8' )  then
            skyCover(i,j)= 'OVC'
         else
            skyCover(i,j)= '   '
         endif
400   continue

c               -------      deal with lacking of data      -------
      do j= 1,staNum
         autoStationType(j)= "UNK"
c        presWeather(j)= "UNK"
         reportType(j)= "SYNOP"
         stationName(j)= "UNK"

         altimeter(j)= badflag
         precip1Hour(j)= badflag
         precip6Hour(j)= badflag
         snowCover(j)= badflag
      enddo

1000  return
      end




      subroutine read_meso_cwb (inpath, maxobs, badflag, ibadflag,
     ~                          i4time_sys, timeObs, rptTp, stnTp, 
     ~                          stname, lats, lons, elev,
     ~                          t, t24max, t24min, td, rh, 
     ~                          pcp1hr, pcp3hr, pcp6hr, pcp24hr, 
     ~                          dd, ff, wgdd, wgff,
     ~                          stnp, mslp, pcc, pc, sr, st,
     ~                          num, istatus)                    
 
c======================================================================
c
c     Routine to read the CWB ASCII Mesonet files.
c     
c======================================================================
 
      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real :: rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real :: pcp24hr(maxobs), dd(maxobs), ff(maxobs), wgdd(maxobs)
      real :: wgff(maxobs), stnp(maxobs), mslp(maxobs), pc(maxobs)
      real :: sr(maxobs), st(maxobs)
      integer :: pcc(maxobs), wmoId(maxobs)

      double precision  timeObs(maxobs)

      logical :: l_parse

c    larger arrays for istart and iend to read data to make processes smooth
      integer, parameter :: num40 = 40,  num70 = 70
      integer   :: istart(num70), iend(num70), hh, flag
 
      character(*)  :: inpath
      character(13) :: cvt_i4time_wfo_fname13, a13time_eat
      character(6)  :: rptTp(maxobs), stnTp(maxobs)
      character(5)  :: stname(maxobs), c5_blank
      character(3)  :: cstn_id, stn_id(maxobs)
      character     :: filename*80, line*320
 
c                      Stuff for the mesonet metadata.
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)
 
      character :: stn_id_master(maxobs)*3, stn_name_master(maxobs)*5
 
c               Get the mesonet metadata (station information).
      call read_tmeso_stntbl (inpath, maxobs, badflag,  
     ~                        stn_id_master, stn_name_master,
     ~                        lat_master, lon_master, elev_master,
     ~                        num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' Error reading mesonet station table'
         return
      endif

c    Fill the output arrays with something, then open the file to read.
 
      istatus=  0
      c5_blank= '     '
      rptTp=    'LDAD'
      stnTp=    'MESO'
      stname=   c5_blank 
      pcc =     ibadflag
      t   =     badflag
      td  =     badflag
      rh  =     badflag
      stnp=     badflag
      dd  =     badflag
      ff  =     badflag
      wgdd=     badflag
      wgff=     badflag
      pc  =     badflag
      sr  =     badflag
      st  =     badflag
      timeObs=  dble( i4time_sys )

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)

      filename= 'Data.CWB.MSO.'
     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)            ! yyyy_mm
     ~           //'-'//a13time_eat(7:8) //'_' //'0000' //'_h.pri'    ! dd
      write(6,*) ' Mesonet file ', filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0
      num_keep= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname), 
     ~         status='old',err=980)

c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 100  flag= 0
 
      read (11,'(a)',end=600,err=990) line

c    Find first dash in time portion (two spaces before last dash in string)
      do i= 1,300
         if (line(i:i) == ':')  exit
      enddo 
 
c                  Parse the string into contiguous characters
      idash= i -9
      istart= 0
      iend=   0

      ivar= 1
      istart(1)= 1

      do i= 1,idash
         if ( i == 1 )  go to 200
         if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then
            iend(ivar)= i-1
         endif

 200     if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
            ivar= ivar +1
            istart(ivar)= i+1
         endif
      enddo

      if ( istart(num40) /= 0 )  go to 100

      ivar= 1
      read (line(istart(ivar):iend(ivar)),'(2i2)',err=399) ihr, imin

      read (a13time_eat(10:11),'(i2)') hh
      if ( ihr /= hh )  go to 100

      ivar= 2
      read (line(istart(ivar):iend(ivar)),*,err=399) cstn_id

      ivar= 3
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rstnp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rstnp
      endif

      ivar= 4
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         slp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) slp 
      endif

      ivar= 5
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         ipcc= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) ipcc
      endif

      ivar= 6
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpc= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rpc
      endif

      ivar= 7
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rt= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rt
      endif

      ivar= 10
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rtd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rtd
      endif

      ivar= 11
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         idir= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) idir
      endif

      ivar= 12
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rspd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rspd
      endif

      ivar= 13
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rwgff= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rwgff
      endif

      ivar= 14
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         iwgdd= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) iwgdd
      endif

      ivar= 15
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpcp= badflag
      elseif ( l_parse(line(istart(ivar):iend(ivar)),'000T') ) then
         rpcp= 0
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rpcp
      endif

      ivar= 21
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rsr= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rsr
      endif

      ivar= 24
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         irh= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) irh
      endif

      ivar= 28
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rst= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rst
      endif

      if ( num == 0 )  go to 400

      do i= 1,num
         if ( cstn_id//'  ' == stn_id(i) ) then 
            num= i
            flag= 1
            go to 500
         else
            cycle
         endif
      enddo
      go to 400

 399  write(6,*) ' read error in station/variable ', num+1, ivar
      write(6,*) ivar, line(istart(ivar):iend(ivar))
      go to 990

c    Have good date/time...store ob.  Adjust/scale variables while storing.
 400  num= num_keep +1                    ! add to count

      if ( num > maxobs ) then
         write(6,*) ' read_local_cwb error for too many obs: ',
     ~              num, maxobs
         istatus= 0
         return
      endif
 
c Match data with metadata for this station, then store the metadata in arrays.
      imatch= 0
      do j= 1,num_master
         if ( cstn_id == stn_id_master(j) ) then
            lats(num)= lat_master(j)
            lons(num)= lon_master(j)
            elev(num)= elev_master(j)
            stn_id(num)= stn_id_master(j) 
            stname(num)= stn_name_master(j) 
            imatch=1
         endif
      enddo 

      if ( imatch == 0 ) then
         write(6,*) ' No station match ', cstn_id
      endif
 
c     stname(num)= stn_id//'  '
 
c                                quality control
 500  if ( rstnp <= 0 ) then
         stnp(num)= badflag
      else
         stnp(num)= rstnp *100.                        ! millibar -> pascal
      endif
 
      if ( slp <= 800. .or. slp > 1100. ) then
         mslp(num)= badflag
      else
         mslp(num)= slp *100.                          ! millibar -> pascal
      endif
 
      pcc(num)= ipcc 

      if ( rpc <= 0 ) then
         pc(num)= badflag
      else
         pc(num)= rpc *100.                            ! millibar -> pascal 
      endif
 
      if ( rt <= -90 ) then
         t(num)= badflag
      else
         if ( rt > 50. )  rt= - (rt - 50.)
         t(num)= rt +273.15                            ! degC -> degK
      endif
 
      if ( rtd <= -90 ) then
         td(num)= badflag
      else
         if ( rtd > 50. ) rtd= - (rtd - 50.)
         td(num)= rtd +273.15                          ! degC -> degK
      endif
 
      if ( idir > 36 .or. idir < 0 ) then
         dd(num)= badflag
      else
         dd(num)= float(idir * 10)                     ! unit : deg
      endif
 
      if ( rspd < 0 ) then
         ff(num)= badflag
      else
         ff(num)= rspd                                 ! unit : m/s 
      endif
 
      if ( rwgff < 0 ) then
         wgff(num)= badflag
      else
         wgff(num)= rwgff                              ! unit : m/s
      endif
 
      if ( iwgdd > 36 .or. iwgdd < 0 ) then
         wgdd(num)= badflag
      else
         wgdd(num)= float(iwgdd * 10)                  ! unit : deg
      endif
 
      if ( rpcp < 0 ) then
         pcp1hr(num)= badflag
      else
         pcp1hr(num)= rpcp *0.001                      ! millimeter -> meter
      endif
 
      if ( rsr < 0 ) then
         sr(num)= badflag
      else
         sr(num)= rsr /1000. /3600.                    ! conv mJ/m/m to watt/m/m
      endif
 
      if ( irh < 0 ) then
         rh(num)= badflag
      else
         rh(num)= float(irh)                           ! unit : %
      endif
 
      if ( rst < 0 ) then
         st(num)= badflag
      else
         st(num)= rst                                  ! degC -> degK
      endif
 
c                          Go back for the next ob.
      if ( flag == 0 )  num_keep= num 
      num= num_keep
      go to 100
 
 600  call mso_t24_pcp (inpath, filename, stname, stn_id, maxobs, 
     ~     badflag, hh, num, t24max, t24min, pcp3hr, pcp6hr, pcp24hr,      
     ~     istatus)

      if ( istatus == 1 ) then
         do i= 1,maxobs
            t24max(i)= t24max(i) +273.15                  ! degC -> degK
            t24min(i)= t24min(i) +273.15                  ! degC -> degK
            pcp3hr(i)= pcp3hr(i) *0.001                   ! millimeter -> meter
            pcp6hr(i)= pcp6hr(i) *0.001                   ! millimeter -> meter
            pcp24hr(i)= pcp24hr(i) *0.001                 ! millimeter -> meter
         enddo
      else
         write(6,*) ' Error estimating mso_t24_pcp '
      endif

c                        Hit end of file...that's it.
      write(6,*) ' Found ', num, ' mesonet stations.'
      istatus= 1
      return
      
 980  write(6,*) ' Warning: could not open mesonet data file ',filename
      num= 0
      istatus= -1
      return

 990  write(6,*) ' ** ERROR reading mesonet data.'
      num= 0
      istatus= -1
      return
      
      end
 
 
 
      subroutine read_tmeso_stntbl (inpath, maxobs, badflag, stn_id,
     ~           stn_name, lat, lon, elev, num, istatus)       
 
c======================================================================
c
c     Routine to read station information for the CWB ASCII Mesonet 
c	data.
c     
c======================================================================
 
      real         :: lat(maxobs), lon(maxobs), elev(maxobs)
      character(3) :: stn_id(maxobs), stn_id_in
      character(5) :: stn_name(maxobs), stn_name_in
      character(*) :: inpath
 
      lat=  badflag
      lon=  badflag
      elev= badflag
      stn_id=   '   '
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'stn-table',status='old',
     ~                                                err=990)
      num= 0

c                 Skip header comments at the top of the file
      do iread= 1,2
         read (13,*,end=550,err=990)
      enddo
 
c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.

 500  read (13,900,end=550,err=990) stn_id_in,stn_name_in,
     ~                              lat_deg,lat_min,lat_sec,alat_sec,       
     ~                              lon_deg,lon_min,lon_sec,alon_sec,
     ~                              elev_m
 900  format (2x,a3,1x,a5,14x,                     ! name
     ~        i2,2x,i2,1x,i2,1x,f3.0,4x,           ! lat
     ~        i3,2x,i2,1x,i2,1x,f3.0,              ! lon
     ~        f12.0)                               ! elevation
 
c         Move station info to arrays for sending to calling routine.
 
      alat= float(lat_deg) +float(lat_min)/60. 
     ~                     +(float(lat_sec) +alat_sec) /3600.
      alon= float(lon_deg) +float(lon_min)/60. 
     ~                     +(float(lon_sec) +alon_sec) /3600.

      num= num +1
      stn_id(num)= stn_id_in
      stn_name(num)= stn_name_in
      lat(num)= alat
      lon(num)= alon
      elev(num)= elev_m
 
c                         Go back for the next ob.
      go to 500
 
c                        Hit end of file...that's it.
 550  write(6,*) ' Found ', num,
     ~           ' mesonet stations in the station table.' 
      istatus= 1
      return
      
 980  write(6,*) ' Warning: could not open mesonet station file ',
     ~           inpath
      istatus= -1
      return

 990  write(6,*) stn_id_in, stn_name_in,
     ~           lat_deg, lat_min, lat_sec, alat_sec,       
     ~           lon_deg, lon_min, lon_sec, alon_sec, elev_m
      write(6,*) ' ** ERROR reading mesonet station table'
      istatus= 0
      return

      end


 
      subroutine mso_t24_pcp (inpath, filename, stname, stn_id, maxobs,
     ~           badflag, ih, num, t24max, t24min, pcp3hr, pcp6hr, 
     ~           pcp24hr, istatus)

      integer, parameter :: num24 = 24,  num40 = 40,  num70 = 70
      character(*) :: stname(maxobs), stn_id(maxobs), inpath
      character(2) :: yy, mm, dd
      character    :: filename*35, fileDummy*35, line*320, stn*3
      logical :: l_parse
      integer :: istart(num70), iend(num70), d(12) 
      integer :: hr, flag
      real :: t24max(maxobs), t24min(maxobs)
      real :: tmax(maxobs,num24), tmin(maxobs,num24)
      real :: pcp3hr(maxobs), pcp6hr(maxobs), pcp24hr(maxobs)
      real :: p1hr(maxobs,-4:num24)
      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      stn=  '   '
      tmax= badflag
      tmin= badflag
      p1hr=    0
      pcp3hr=  0
      pcp6hr=  0
      pcp24hr= 0
      maxrcd= 3000
      istatus= 0

      read (filename(14:17),'(i4)') iy
      read (filename(19:20),'(i2)') im
      read (filename(22:23),'(i2)') id

c                 open two files to read data of 24 hours
      do l= 0,1
         id= id -l

         if ( id < 1 ) then
            im= im -1
  
            if ( im < 1 ) then
               im= 12
               iy= iy -1
            endif

            id= d(im)
         endif

         iy= iy -2000
         call i2a ( iy, yy )
         call i2a ( im, mm )
         call i2a ( id, dd )
         fileDummy= 'Data.CWB.MSO.' //'20' //yy //'-' //mm //'-' //dd
     ~                                          //'_' //'0000_h.pri'     
         iy= iy +2000

         call s_len ( inpath, len_inpath )
         call s_len ( fileDummy, len_fname )

         select case ( l ) 
         case ( 0 )
            rewind (11)
         case ( 1 )
            open (11,file=inpath(1:len_inpath)//fileDummy(1:len_fname),
     ~               status='old',err=980)
         end select

         do 200 k= 1,maxrcd
            read (11,'(a)',iostat=istat) line

            if ( istat == -1 .or. istat == -2 )  exit
            if ( istat > 0 )  go to 990
               
c    Find first dash in time portion (two spaces before last dash in string)
            do i= 1,300
               if (line(i:i) == ':')  exit
            enddo

            idash= i -9

c                 Parse the string into contiguous characters
            istart= 0
            iend=   0

            ivar= 1
            istart(1)= 1

            do i= 1,idash
               if ( i == 1 )  go to 100

               if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then      
                  iend(ivar)= i-1
               endif

 100           if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
                  ivar= ivar +1
                  istart(ivar)= i+1
               endif
            enddo

c                              avoid daily data
            if ( istart(num40) /= 0 )  cycle

            ivar= 1
            read (line(istart(ivar):iend(ivar)),'(2i2)',err=199) hr, mn
            if ( l == 0 .and. hr >  ih )  cycle  
            if ( l == 1 .and. hr <= ih )  cycle  
            if ( hr == 0 )  hr= 24
           
            ivar= 2
            read (line(istart(ivar):iend(ivar)),*,err=200) stn

            do i= 1,num
            if ( stn_id(i) == stn//'  ' ) then
               ivar= 8
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmax(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmax(i,hr)
               endif

               ivar= 9
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmin(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmin(i,hr)    
               endif

               ivar= 15
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 p1hr(i,hr)= badflag
               elseif(l_parse(line(istart(ivar):iend(ivar)),'000T'))then
                 p1hr(i,hr)= 0
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)p1hr(i,hr)
               endif
            endif
            enddo

            cycle

 199        write(6,*) ' read error in station/variable ', j, ivar
            write(6,*) ivar,line(istart(ivar):iend(ivar))
 200     enddo

         if ( ih == 23 )  exit 
 500  enddo 

      do i= 1,num
         t24max(i)= tmax(i,1)
         t24min(i)= tmin(i,1)

         do j= 2,num24
            flag= 0

            if ( tmax(i,j) > 50. .or. tmax(i,j) < -90. ) then
               tmax(i,j)= badflag
               exit
            elseif ( tmin(i,j) > 50. .or. tmin(i,j) < -90. ) then
               tmin(i,j)= badflag
               exit
            endif

            if ( t24max(i) < tmax(i,j) )  t24max(i)= tmax(i,j)
            if ( t24min(i) > tmin(i,j) )  t24min(i)= tmin(i,j)

            flag= 1
         enddo

c  once there is any data missing in tmax or tmin, flag= 0 from prior do loop
         if ( flag /= 1 ) then
            t24max(i)= badflag
            t24min(i)= badflag
            write(6,*) 'too few data to obtain Tmax/Tmin for ',
     ~                 stname(i), ' mesonet station ', j
            cycle
         endif
 900  enddo

c                         calculate accumulated rain gauge
      do i= 1,num
         p1hr(i, 0)= p1hr(i,24)
         p1hr(i,-1)= p1hr(i,23)
         p1hr(i,-2)= p1hr(i,22)
         p1hr(i,-3)= p1hr(i,21)
         p1hr(i,-4)= p1hr(i,20)
      enddo

      do i= 1,num
      do j= ih,ih-2,-1
         if ( p1hr(i,j) == badflag ) then
            pcp3hr(i)= badflag
            write(6,*) 'too few data to estimate pcp3hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp3hr(i)= pcp3hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= ih,ih-5,-1
         if ( p1hr(i,j) == badflag ) then
            pcp6hr(i)= badflag
            write(6,*) 'too few data to estimate pcp6hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp6hr(i)= pcp6hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= 1,num24
         if ( p1hr(i,j) == badflag ) then
            pcp24hr(i)= badflag
            write(6,*) 'too few data to estimate pcp24hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp24hr(i)= pcp24hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      istatus= 1
      return

 980  write(6,*) ' Warning: could not open mesonet data file ',
     ~           inpath(1:len_inpath)//fileDummy(1:len_fname)
      istatus= -1
      return
        
 990  write(6,*) ' ** ERROR reading mesonet data.'
      num= 0
      istatus= -1
      return

      end
        


      subroutine  i2a (ii,aa)

      character(2) :: aa
      integer      :: ii

      if ( ii < 10 ) then
         write (aa,'(a1,i1)') '0', ii
      else
         write (aa,'(i2)') ii
      endif

      return
      end
