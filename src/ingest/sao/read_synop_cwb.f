      subroutine read_synop_cwb (
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
      character*6    autoStationType(recNum)
      character*25   presWeather(recNum)
      character*6    reportType(recNum)
      character*8    skyCover(maxSkyCover,recNum)
      character*5    stationName(recNum)

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

      character*3   reportFlag(recNum)
      character*2   yy(recNum), mo(recNum), dd(recNum)
      character*2   hh(recNum), mn(recNum)
      character*10  time(recNum)
      character*9   a10_to_a9

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
         if(dummy .gt. 0)then
           if ( skyCover(i,j) .eq. skyCover(1,j)  .and.
     ~        skyLayerBase(i,j) .ge. lowestCloudHeight(dummy) )  then
            if ( dummy .ne. 9  .and.
     ~           skyLayerBase(i,j) .gt. lowestCloudHeight(dummy+1) )
     ~         go to 250
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
         presWeather(j)= "UNK"
         reportType(j)= "SYNOP"
         stationName(j)= "UNK"

         altimeter(j)= badflag
         precip1Hour(j)= badflag
         precip6Hour(j)= badflag
         snowCover(j)= badflag
      enddo

1000  return
      end
