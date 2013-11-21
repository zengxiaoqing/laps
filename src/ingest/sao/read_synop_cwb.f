      subroutine read_synop_cwb ( filename, maxSkyCvr, maxobs,
     &                          i4time_sys, path_to_local,            
     &                          altm, stnTp, td, tdTths, elev,
     &                          lats, lons, t24max, t24min,
     &                          pcp1hr, pcp24hr, pcp3hr, pcp6hr, 
     &                          prsWth, p, pc, pcc, rptTp, rh,
     &                          mslp, skyCvr, skyLyrBs, snowCvr, sr, st,       
     &                          stname, tTths, t, timeObs, vis,
     &                          dd, wgdd, ff, wgff, wmoId, badflag,
     &                          num, istatusSynop )

      integer, parameter :: maxSynop = 150
      integer, parameter :: maxMso =    60

      character*(*)  filename
      character(*)   path_to_local
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
     &     altm(1:nq), stnTp(1:nq), td(1:nq), tdTths(1:nq), elev(1:nq),
     &     lats(1:nq), lons(1:nq), t24max(1:nq), t24min(1:nq),
     &     pcp1hr(1:nq), pcp24hr(1:nq), pcp3hr(1:nq), pcp6hr(1:nq), 
     &     prsWth(1:nq), pc(1:nq), pcc(1:nq), rptTp(1:nq), mslp(1:nq),
     &     skyCvr(1:maxSkyCvr,1:nq), skyLyrBs(1:maxSkyCvr,1:nq),
     &     snowCvr(1:nq), stname(1:nq), tTths(1:nq), t(1:nq),
     &     timeObs(1:nq), vis(1:nq), dd(1:nq), wgff(1:nq), ff(1:nq),
     &     wmoId(1:nq), badflag, numSynop, istatusSynop )

      do i= 1,numSynop
         write(stnNo(i),'(i5)') wmoId(i)      
         stname(i)(1:5)=stnNo(i)(1:5)
      enddo

      np= numSynop +1
      nq= numSynop +maxMso

      call s_len ( path_to_local, len_inpath )
      path_to_local= path_to_local(1:len_inpath)//'mso/'

      call read_meso_cwb ( path_to_local, maxMso, badflag, ibadflag, 
     &                     i4time_sys, timeObsMso, rptTpMso, stnTpMso,
     &                     stnNoMso, latsMso, lonsMso, elevMso,
     &                     tMso, t24maxMso, t24minMso, tdMso, rhMso, 
     &                     pcp1hrMso, pcp3hrMso, pcp6hrMso, pcp24hrMso,
     &                     ddMso, ffMso, wgddMso, wgffMso, pMso, 
     &                     mslpMso, pccMso, pcMso, srMso, stMso, 
     &                     numMso, istatusMso )

c                    combine synop data and mesonet data 
      do i= 1,numMso ! added by Shuyuan20100707 for checking meso RH values
         if(rhMso(i)>101) then
             rhMso(i)=  badflag       
         endif
      enddo
c

      do i= 1,numSynop ! added by Shuyuan20100707 for checking synop RH data 
         if(rh(i)>101) then
             rh(i)=  badflag       
         endif
      enddo

      
      k= numSynop
      do i= 1,numMso
         flag= 0

         do j= 1,numSynop
            if ( stnNoMso(i) == stnNo(j) ) then
               rh(j)=     rhMso(i)
               dd(j)=     ddMso(i)
               ff(j)=     ffMso(i)
               wgdd(j)=   wgddMso(i)
               wgff(j)=   wgffMso(i)
               p(j)=      pMso(i)
               sr(j)=     srMso(i)
               st(j)=     stMso(i)
c
c modified by min-ken.hsieh
c must assign timeObsMso here
c or timeObs will be synop obs time
c Also, assign mso t/td to synop t/td and tTths/tdTths
c or we will get t/td carried with synop which only update hourly
c
	       t(j)=	  tMso(i)
	       tTths(j)=  tMso(i)
	       td(j)=	  tdMso(i)
	       tdTths(j)= tdMso(i)
               timeObs(j)=timeObsMso(i)
	       pcp1hr(j)= pcp1hrMso(i)
	       pcp3hr(j)= pcp3hrMso(i)
	       pcp6hr(j)= pcp6hrMso(i)
	       pcp24hr(j)= pcp24hrMso(i)
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
            stname(k)=  stnNoMso(i)
            wmoId(k)=   ibadflag
            td(k)=      tdMso(i)
            tdTths(k)=  tdMso(i)
            elev(k)=    elevMso(i)
            lats(k)=    latsMso(i)
            lons(k)=    lonsMso(i)
            t24max(k)=  t24maxMso(i)
            t24min(k)=  t24minMso(i)
            rh(k)=      rhMso(i)

            pcp1hr(k)=  pcp1hrMso(i)
            pcp24hr(k)= pcp24hrMso(i)
            pcp3hr(k)=  pcp3hrMso(i)
            pcp6hr(k)=  pcp6hrMso(i)
            p(k)=       pMso(i)
c         add if shuyuan  20100707 
            if(pcMso(i)>0.5  .and. pcMso(i)<1000)  then
              pc(k)=  pcMso(i)
            else
            pc(k)=badflag
            endif 
            if(pccMso(i)>-10000  .and. pccMso(i)<10000)  then   
              pcc(k)=     pccMso(i)
            else
              pcc(k)=badflag
            endif
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

      num = k

      end



      subroutine read_synop_cwb_sub (
     &     filename, maxSkyCover, recNum, altimeter,
     &     autoStationType, dewpoint, dpFromTenths, elevation,
     &     latitude, longitude, maxTemp24Hour, minTemp24Hour,
     &     precip1Hour, precip24Hour, precip3Hour, precip6Hour,
     &     presWeather, pressChange3Hour, pressChangeChar,
     &     reportType, seaLevelPress, skyCover, skyLayerBase,
     &     snowCover, stationName, tempFromTenths, temperature,
     &     timeObs, visibility, windDir, windGust, windSpeed, wmoId,
     &     badflag, staNum, istatus )

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
      character(8)   skyCoverDmy(recNum)

      integer  windQua(recNum), seaLevelPressQua(recNum)
      integer  temperatureQua(recNum), dewpointQua(recNum)
      integer  pressChange3HourQua(recNum)
      integer  dupliStation(9), staNum, dummy

c     real  tempDewDiff(recNum), lowestCloudHeight(0:10)
c     real  skyLayerBaseDmy(recNum)
      real  skyLayerBaseDmy(recNum), tempDewDiff(recNum)
      real  lowestCloudHeight(0:10)

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
     ~                   skyCoverDmy(j), skyLayerBaseDmy(j)
         read (1,40,end=9,err=9) tempDewDiff(j), dewpointQua(j),       
     ~                   pressChangeChar(j), pressChange3Hour(j),
     ~                   pressChange3HourQua(j), precip3Hour(j),
     ~                   maxTemp24Hour(j), minTemp24Hour(j), windGust(j) 
         read (1,50,end=9,err=9) skyCover(1,j), skyLayerBase(1,j),
     ~                   skyCover(2,j), skyLayerBase(2,j), 
     ~                   skyCover(3,j), skyLayerBase(3,j), 
     ~                   precip24Hour(j)
         read (1,60,end=9,err=9) skyCover(4,j), skyLayerBase(4,j)

         if ( reportFlag(j) /= '*31' )  then
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
     ~               skyCoverDmy(j), skyLayerBaseDmy(j)
         write (6,*) dewpoint(j), dewpointQua(j),
     ~               pressChangeChar(j), pressChange3Hour(j),
     ~               pressChange3HourQua(j), precip3Hour(j),
     ~               maxTemp24Hour(j), minTemp24Hour(j), windGust(j) 
         write (6,*) skyCover(1,j), skyLayerBase(1,j),
     ~               skyCover(2,j), skyLayerBase(2,j), 
     ~               skyCover(3,j), skyLayerBase(3,j), precip24Hour(j)
         write (6,*) skyCover(4,j), skyLayerBase(4,j)

10       staNum= staNum +1
      enddo

20    format ( a3, i5, f4.0, 2f5.2, 2x, 5a2 )
30    format ( 2x, 2f3.0, i1, f3.0, a2, 3x, f5.1, i1, f4.1, i1, a2, 2x,
     ~         f3.0 )
40    format ( f4.1, i1, 1x, i2, f4.1, i1, 3(1x, f4.1), 8x, f3.0 )
50    format ( 3(a2, 2x, f3.0), 8x, f4.1 )
!HJ: W>=D+3. 10/14/2013
60    format ( a2, 2x, f3.0 )

c     --- eliminate duplicate data coming from international broadcast ---
99    do 100 k= 1,9
      do 100 j= 1,staNum
         if ( wmoId(j) == dupliStation(k) )  then
            staNum= staNum -1
            do i= j,staNum
               reportFlag(i)= reportFlag(i+1)
	       wmoId(i)     = wmoId(i+1)
	       elevation(i) = elevation(i+1)
	       latitude(i)  = latitude(i+1)
	       longitude(i) = longitude(i+1)
	       yy(i)= yy(i+1)
	       mo(i)= mo(i+1) 
	       dd(i)= dd(i+1)
	       hh(i)= hh(i+1)
	       mn(i)= mn(i+1)
	       windDir(i)            = windDir(i+1)
	       windSpeed(i)          = windSpeed(i+1)
	       windQua(i)            = windQua(i+1)
	       visibility(i)         = visibility(i+1)
	       presWeather(i)        = presWeather(i+1)
	       seaLevelPress(i)      = seaLevelPress(i+1)
	       seaLevelPressQua(i)   = seaLevelPressQua(i+1)
	       temperature(i)        = temperature(i+1)
	       temperatureQua(i)     = temperatureQua(i+1)
               skyCoverDmy(i)        = skyCoverDmy(i+1)
               skyLayerBaseDmy(i)    = skyLayerBaseDmy(i+1)
	       tempDewDiff(i)        = tempDewDiff(i+1)
	       dewpointQua(i)        = dewpointQua(i+1)
	       pressChangeChar(i)    = pressChangeChar(i+1)
	       pressChange3Hour(i)   = pressChange3Hour(i+1)
	       pressChange3HourQua(i)= pressChange3HourQua(i+1)
	       precip3Hour(i)        = precip3Hour(i+1)
	       maxTemp24Hour(i)      = maxTemp24Hour(i+1)
	       minTemp24Hour(i)      = minTemp24Hour(i+1)
	       windGust(i)           = windGust(i+1)
	       precip24Hour(i)       = precip24Hour(i+1)
               do l= 1,maxSkyCover
	          skyCover(l,i)=     skyCover(l,i+1)
		  skyLayerBase(l,i)= skyLayerBase(l,i+1)
	       enddo	
            enddo
         endif
100   continue
      write(*,*) 'SYNOP staNum=', staNum

c      ----------       examine data quality and change units       ---------
      do j= 1,staNum
         if ( windQua(j) /= 1 )  then
            windDir(j)= badflag
            windSpeed(j)= badflag
         endif

         if ( windGust(j) == -99. )  windGust(j)= badflag
         if ( elevation(j) == -999. )  elevation(j)= badflag
         if ( pressChangeChar(j)==-9 ) pressChangeChar(j)=int(badflag)
         if ( presWeather(j) == '-9' )  presWeather(j)= 'UNK'

         if ( seaLevelPressQua(j) == 1 )  then
               seaLevelPress(j)= seaLevelPress(j) *100.   ! millibar -> pascal
            else
               seaLevelPress(j)= badflag
         endif

         if ( pressChange3HourQua(j) == 1 )  then
            pressChange3Hour(j)= pressChange3Hour(j) *100. ! millibar -> pascal
          else
            pressChange3Hour(j)= badflag
            pressChangeChar(j)= int(badflag)
         endif

         if ( temperatureQua(j) == 1 )  then
               temperature(j)= temperature(j) +273.15           ! degC -> degK
            else
               temperature(j)= badflag
         endif

         tempFromTenths(j)= temperature(j)

         if ( dewpointQua(j) == 1 )  then
               dewpoint(j)= temperature(j) -tempDewDiff(j)        ! unit: degK
            else
               dewpoint(j)= badflag
         endif

         dpFromTenths(j)= dewpoint(j)

         if ( maxTemp24Hour(j) == -99.9 )  then
               maxTemp24Hour(j)= badflag
            else
               maxTemp24Hour(j)= maxTemp24Hour(j) +273.15       ! degC -> degK
         endif

         if ( minTemp24Hour(j) == -99.9 )  then
               minTemp24Hour(j)= badflag
            else
               minTemp24Hour(j)= minTemp24Hour(j) +273.15       ! degC -> degK
         endif

         if ( precip24Hour(j) == -99.9 )  then
               precip24Hour(j)= badflag
            else
               precip24Hour(j)= precip24Hour(j) *0.001   ! millimeter -> meter
         endif

         if ( precip3Hour(j) == -99.9 )  then
               precip3Hour(j)= badflag
            else
               precip3Hour(j)= precip3Hour(j) *0.001     ! millimeter -> meter
         endif

         if ( yy(j)(1:1) == ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( mo(j)(1:1) == ' ' )  mo(j)= '0'//mo(j)(2:2)
         if ( dd(j)(1:1) == ' ' )  dd(j)= '0'//dd(j)(2:2)
         if ( hh(j)(1:1) == ' ' )  hh(j)= '0'//hh(j)(2:2)
         if ( mn(j)(1:1) == ' ' )  mn(j)= '0'//mn(j)(2:2)
         time(j)= yy(j)//mo(j)//dd(j)//hh(j)//mn(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeObs(j)= dble( i4time )                       ! seconds since 1960
      enddo

c    -------    transform code figure into visibility ( unit: m )  -------
      do j= 1,staNum
         if     ( visibility(j) ==  0. )  then
               visibility(j)= 50.
         elseif ( visibility(j) >  0.  .and. 
     ~            visibility(j) < 51. )  then
               visibility(j)= visibility(j) *100.
         elseif ( visibility(j) > 55.  .and. 
     ~            visibility(j) < 81. )  then
               visibility(j)= ( visibility(j) -50. ) *1000.
         elseif ( visibility(j) > 80.  .and. 
     ~            visibility(j) < 90. )  then
               visibility(j)= ( visibility(j) -74. ) *5000.
         elseif ( visibility(j) == 90. )  then
               visibility(j)= 25.
         elseif ( visibility(j) == 91. )  then
               visibility(j)= 50.
         elseif ( visibility(j) == 92. )  then
               visibility(j)= 200.
         elseif ( visibility(j) == 93. )  then
               visibility(j)= 500.
         elseif ( visibility(j) == 94. )  then
               visibility(j)= 1000.
         elseif ( visibility(j) == 95. )  then
               visibility(j)= 2000.
         elseif ( visibility(j) == 96. )  then
               visibility(j)= 4000.
         elseif ( visibility(j) == 97. )  then
               visibility(j)= 10000.
         elseif ( visibility(j) == 98. )  then
               visibility(j)= 20000.
         elseif ( visibility(j) == 99. )  then
               visibility(j)= 50000.
         else
               visibility(j)= badflag
         endif
      enddo

c ---- transform code figure into base of cloud layer indicated (unit: m) ----
      do 200 i= 1,maxSkyCover
      do 200 j= 1,staNum
         if     ( skyLayerBase(i,j) ==  0. )  then
               skyLayerBase(i,j)= 15.
         elseif ( skyLayerBase(i,j) >  0.  .and.  
     ~            skyLayerBase(i,j) < 51. )  then
               skyLayerBase(i,j)= skyLayerBase(i,j) *30.
         elseif ( skyLayerBase(i,j) > 55.  .and.  
     ~            skyLayerBase(i,j) < 81. )  then
               skyLayerBase(i,j)= ( skyLayerBase(i,j) -50. ) *300.
         elseif ( skyLayerBase(i,j) > 80.  .and. 
     ~            skyLayerBase(i,j) < 90. )  then
               skyLayerBase(i,j)= ( skyLayerBase(i,j) -74. ) *1500.
         elseif ( skyLayerBase(i,j) == 90. )  then
               skyLayerBase(i,j)= 25.
         elseif ( skyLayerBase(i,j) == 91. )  then
               skyLayerBase(i,j)= 75.
         elseif ( skyLayerBase(i,j) == 92. )  then
               skyLayerBase(i,j)= 150.
         elseif ( skyLayerBase(i,j) == 93. )  then
               skyLayerBase(i,j)= 250.
         elseif ( skyLayerBase(i,j) == 94. )  then
               skyLayerBase(i,j)= 450.
         elseif ( skyLayerBase(i,j) == 95. )  then
               skyLayerBase(i,j)= 800.
         elseif ( skyLayerBase(i,j) == 96. )  then
               skyLayerBase(i,j)= 1250.
         elseif ( skyLayerBase(i,j) == 97. )  then
               skyLayerBase(i,j)= 1750.
         elseif ( skyLayerBase(i,j) == 98. )  then
               skyLayerBase(i,j)= 2250.
         elseif ( skyLayerBase(i,j) == 99. )  then
               skyLayerBase(i,j)= 3000.
         else
               skyLayerBase(i,j)= badflag
         endif
200   continue

      do 300 j= 1,staNum
         if (skyCover(1,j) == '-9' .or. skyLayerBase(1,j) == -9.)  cycle       

c          ----- eliminate duplicate skyCovers and skyLayerBases -----
         dummy= int( skyLayerBaseDmy(j) )
         skyLayerBaseDf= abs(skyLayerBase(1,j)-lowestCloudHeight(dummy))     
         if ( dummy >= 0 )  then
            if ( skyCover(1,j) <= skyCoverDmy(j)  .and.
     ~           (skyLayerBase(1,j) >= lowestCloudHeight(dummy) .or.
     ~            skyLayerBaseDf < 30.) )  then
               if ( dummy /= 9  .and.
     ~              skyLayerBase(1,j) > lowestCloudHeight(dummy+1) )
     ~            go to 250
               skyCoverDmy(j)= '   '
	       skyLayerBaseDmy(j)= badflag
            endif
         endif

c        ----- eliminate unreasonable skyCovers and skyLayerBases -----
250      if ( skyCover(1,j)==' 8' .and. skyCoverDmy(j)==' 8' )  then
            skyCoverDmy(j)= '   '
            skyLayerBaseDmy(j)= badflag
         endif
300   continue

c  -----  transform code figure into base of lowest cloud ( unit: m )  -----
      do j= 1,staNum
	 if     ( skyLayerBaseDmy(j) == 0. )  then
       	       skyLayerBaseDmy(j)= 25.
         elseif ( skyLayerBaseDmy(j) == 1. )  then
               skyLayerBaseDmy(j)= 75.
	 elseif ( skyLayerBaseDmy(j) == 2. )  then
	       skyLayerBaseDmy(j)= 150.
         elseif ( skyLayerBaseDmy(j) == 3. )  then
	       skyLayerBaseDmy(j)= 250.
	 elseif ( skyLayerBaseDmy(j) == 4. )  then
	       skyLayerBaseDmy(j)= 450.
	 elseif ( skyLayerBaseDmy(j) == 5. )  then
	       skyLayerBaseDmy(j)= 800.
	 elseif ( skyLayerBaseDmy(j) == 6. )  then
	       skyLayerBaseDmy(j)= 1250.
	 elseif ( skyLayerBaseDmy(j) == 7. )  then
	       skyLayerBaseDmy(j)= 1750.
	 elseif ( skyLayerBaseDmy(j) == 8. )  then
	       skyLayerBaseDmy(j)= 2250.
	 elseif ( skyLayerBaseDmy(j) == 9. .and. skyCoverDmy(j) /= ' 0'
     ~             .and. skyCoverDmy(j) /= '-9' )  then
	       skyLayerBaseDmy(j)= 3000.
	 else
	       skyLayerBaseDmy(j)= badflag
  	 endif
      enddo
 
c   assign the lowest cloud data to the first array when the latter is missing 
      do j= 1,staNum
	 if ( skyLayerBaseDmy(j) /= badflag .and. 
     ~        skyLayerBase(1,j) == badflag )  then
            skyLayerBase(1,j)= skyLayerBaseDmy(j)
            skyCover(1,j)=     skyCoverDmy(j)
         endif
      enddo
	    
c        --- transform code figure of cloud cover into metar format ---
      do 400 i= 1,maxSkyCover
      do 400 j= 1,staNum
         if     ( skyCover(i,j) == ' 0' )  then
            skyCover(i,j)= 'SKC'
	    skyLayerBase(i,j)= 22500.
         elseif ( skyCover(i,j) == ' 1'  .or.
     ~            skyCover(i,j) == ' 2' )  then
            skyCover(i,j)= 'FEW'
         elseif ( skyCover(i,j) == ' 3'  .or.
     ~            skyCover(i,j) == ' 4' )  then
            skyCover(i,j)= 'SCT'
         elseif ( skyCover(i,j) == ' 5'  .or.
     ~            skyCover(i,j) == ' 6'  .or.
     ~            skyCover(i,j) == ' 7' )  then
            skyCover(i,j)= 'BKN'
         elseif ( skyCover(i,j) == ' 8' )  then
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
      integer   :: istart(num70), iend(num70), hhmm, hh, flag
 
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
c
c modified by min-ken,hsieh
c filename include hh information
c to read *_m.pri (to get each 15 min data) instead of *_h.pri
c because of the data format is different
c modified each field parsing below
c

      filename= 'Data.CWB.MSO.'
     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)            ! yyyy_mm
     ~           //'-'//a13time_eat(7:11) //'00' //'_m.pri'    ! dd
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

c
c modified by min-ken,hsieh
c because m.pri file data format is different from h.pri
c each column means different data
c
      ivar= 1
      read (line(istart(ivar):iend(ivar)),'(i4)',err=399) ihrmin

      read (a13time_eat(10:13),'(i4)') hhmm
      read (a13time_eat(10:11),'(i2)') hh
      if ( ihrmin /= hhmm )  go to 100

      ivar= 2
      read (line(istart(ivar):iend(ivar)),*,err=399) cstn_id

      ivar= 9
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rstnp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rstnp
      endif

      ivar= 13
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         slp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) slp 
      endif

      ivar= 10
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rt= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rt
      endif

      ivar= 11
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rtd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rtd
      endif

      ivar= 4
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         idir= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) idir
      endif

      ivar= 3
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rspd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rspd
      endif

      ivar= 12
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpcp= badflag
      elseif ( l_parse(line(istart(ivar):iend(ivar)),'000T') ) then
         rpcp= 0
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rpcp
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
 
      ! Shuyuan & Yuanfu: We found rpc and ipcc are un-initialized
      ! and are not read from the meso file. We remove the pcc and pc assignment

c      pcc(num)= ipcc 

c      if ( rpc <= 0 ) then
c         pc(num)= badflag
c      else
c         pc(num)= rpc *100.                            ! millibar -> pascal 
c      endif
 
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
 
      ! Shuyuan & Yuanfu: rwgff/wgdd are never initialized (gust wind)
c      if ( rwgff < 0 ) then
c         wgff(num)= badflag
c      else
c         wgff(num)= rwgff                              ! unit : m/s
c      endif
 
c      if ( iwgdd > 36 .or. iwgdd < 0 ) then
c         wgdd(num)= badflag
c      else
c         wgdd(num)= float(iwgdd * 10)                  ! unit : deg
c      endif
 
      if ( rpcp < 0 ) then
         pcp1hr(num)= badflag
      else
         pcp1hr(num)= rpcp *0.001                      ! millimeter -> meter
      endif
 
      ! Shuyuan & Yuanfu: rsr is never initialized (gust wind)
      if ( rsr < 0 ) then
         sr(num)= badflag
      else
         sr(num)= rsr /1000. /3600.                    ! conv mJ/m/m to watt/m/m
      endif
 
      ! Shuyuan & Yuanfu: We found another un-initialized variable - irh
      ! We have to remove the following for now:
c      if ( irh < 0 ) then
c         rh(num)= badflag
c      else
c         rh(num)= float(irh)                           ! unit : %
c      endif
 
      ! Shuyuan & Yuanfu: rst is never initialized variable
c      if ( rst < 0 ) then
c         st(num)= badflag
c      else
c         st(num)= rst                                  ! degC -> degK
c      endif
 
c                          Go back for the next ob.
      if ( flag == 0 )  num_keep= num 
      num= num_keep
      go to 100
 
 600  call mso_t24_pcp (inpath, filename, stname, stn_id, maxobs, 
     ~     badflag, hh, num, t24max, t24min, pcp3hr, pcp6hr, pcp24hr,      
     ~     istatus)

      if ( istatus == 1 ) then
c
c        modified by min-ken hsieh
c        some stn may not have enough info to calculate t24 and pcp accum.
c        and mso_t24_pcp will return badsfc here
c        so let's check return values before unit conversion
c
         do i= 1,maxobs
	    if (t24max(i).ne.badflag) then
               t24max(i)= t24max(i) +273.15                  ! degC -> degK
	    endif
	    if (t24min(i).ne.badflag) then
               t24min(i)= t24min(i) +273.15                  ! degC -> degK
	    endif
	    if (pcp3hr(i).ne.badflag) then
               pcp3hr(i)= pcp3hr(i) *0.001                   ! millimeter -> meter
	    endif
	    if (pcp6hr(i).ne.badflag) then
               pcp6hr(i)= pcp6hr(i) *0.001                   ! millimeter -> meter
	    endif
	    if (pcp24hr(i).ne.badflag) then
               pcp24hr(i)= pcp24hr(i) *0.001                 ! millimeter -> meter
	    endif
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
c
c modified by min-ken hsieh
c we read _m file to get other data in subroutine read_mso_cwb
c but here we need to get pcp data from _h file
c

c           rewind (11)
            open (11,file=inpath(1:len_inpath)//fileDummy(1:len_fname),
     ~               status='old',err=980)
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
