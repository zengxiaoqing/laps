cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
      subroutine read_ship_cwb ( filename, maxSkyCover, recNum, 
     +     dataPlatformType, dewpoint, elevation, equivWindSpeed10m,
     +     latitude, longitude, precip1Hour, precip24Hour, precip3Hour,        
     +     precip6Hour, presWeather, pressChange3Hour, pressChangeChar,
     +     seaLevelPress, seaSurfaceTemp, skyCover, skyLayerBase, 
     +     stationName, temperature, timeObs, visibility, 
     +     wetBulbTemperature, windDir, windGust, windSpeed, wmoId,
     +     badflag, n, istatus )
 
      integer  recNum

      character*(*)  filename
      character*25   presWeather(recNum)
      character*8    skyCover(maxSkyCover,recNum), stationName(recNum)

      integer  dataPlatformType(recNum), pressChangeChar(recNum)
      integer  wmoId(recNum)

      real  dewpoint(recNum), elevation(recNum)
      real  equivWindSpeed10m(recNum), latitude(recNum)
      real  longitude(recNum), precip1Hour(recNum)
      real  precip24Hour(recNum), precip3Hour(recNum)
      real  precip6Hour(recNum), pressChange3Hour(recNum)
      real  seaLevelPress(recNum), seaSurfaceTemp(recNum)
      real  skyLayerBase(maxSkyCover,recNum), temperature(recNum)
      real  visibility(recNum), wetBulbTemperature(recNum)
      real  windDir(recNum), windGust(recNum), windSpeed(recNum)

      double precision  timeObs(recNum)

      character*3   reportFlag(recNum)
      character*2   yy(recNum), mo(recNum), dd(recNum)
      character*2   hh(recNum), mn(recNum)
      character*10  time(recNum)
      character*9   a10_to_a9

      integer  windQua(recNum), temperatureQua(recNum)
      integer  seaLevelPressQua(recNum), seaSurfaceTempQua(recNum)
      integer  dewpointQua(recNum), pressChange3HourQua(recNum)
      integer  wetBulbTemperatureQua(recNum)

      real  relaHumility(recNum), tempDewDiff(recNum)

      istatus = 0
      n= 0
 
      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do j= 1,recNum
         read ( 1, 10, end=99, err=999 ) reportFlag(j),
     ~               stationName(j), latitude(j), longitude(j),
     ~               yy(j), mo(j), dd(j), hh(j), mn(j)
         read (1,20) windDir(j), windSpeed(j), windQua(j),
     ~               visibility(j), presWeather(j),
     ~               seaLevelPress(j), seaLevelPressQua(j),
     ~               temperature(j), temperatureQua(j),
     ~               skyCover(1,j), skyLayerBase(1,j)
         read (1,30) tempDewDiff(j), dewpointQua(j), pressChangeChar(j),       
     ~               pressChange3Hour(j), pressChange3HourQua(j),
     ~               precip3Hour(j), wetBulbTemperature(j),
     ~               wetBulbTemperatureQua(j), windGust(j)
         read (1,40) skyCover(2,j), skyLayerBase(2,j), precip24Hour(j)
         read (1,50) seaSurfaceTemp(j), seaSurfaceTempQua(j)

         if ( reportFlag(j) .ne. '*34' )  then
            write (6,*) 'read ship data heading error'
            go to 1000
         endif

         n= n+1
      enddo

!HJ: W>=D+3. 10/14/2013
10    format ( a3, a5, 4x, 2f5.2, 2x, 5a2 )
20    format ( 2x, 2f3.0, i1, f23.0, a2, 3x, f5.1, i1, f4.1, i1, a2, 2x,
     ~         f3.0 )
30    format ( f4.1, i1, 1x, i2, f4.1, i1, 1x, f4.1, 10x, f4.1, i1, 3x,
     ~         f3.0 )
40    format ( a2, 2x, f3.0, 14x, f4.1 )
50    format ( 16x, f4.1, i1 )

c      ----------       examing data quality and changing units       ---------
99    do j= 1,n
         if ( windQua(j) .eq. 9 )  then
            windDir(j)= badflag
            windSpeed(j)= badflag
         endif

         if ( windGust(j) .eq. -99. )  windGust(j)= badflag
         if ( elevation(j) .eq. -999. )  elevation(j)= badflag
         if ( precip3Hour(j) .eq. -99.9 )  precip3Hour(j)= badflag
         if ( pressChangeChar(j).eq.-9 ) pressChangeChar(j)=int(badflag)
         if ( presWeather(j) .eq. '-9' )  presWeather(j)= 'UNK'

         if ( seaLevelPressQua(j) .eq. 1 )  then
               seaLevelPress(j)= seaLevelPress(j) *100.   ! millibar -> pascal
            else
               seaLevelPress(j)= badflag
         endif
  
         if ( pressChange3HourQua(j) .eq. 1 )  then       ! millibar -> pascal
               pressChange3Hour(j)= pressChange3Hour(j) *100.
            else
               pressChange3Hour(j)= badflag
         endif

         if ( seaSurfaceTempQua(j) .eq. 1 )  then
               seaSurfaceTemp(j)= seaSurfaceTemp(j) +273.15     ! degC -> degK
            else
               seaSurfaceTemp(j)= badflag
         endif

         if ( temperatureQua(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15           ! degC -> degK
            else
               temperature(j)= badflag
         endif

         if ( dewpointQua(j) .eq. 1 )  then
               dewpoint(j)= temperature(j) -tempDewDiff(j)        ! unit: degK
            else
               dewpoint(j)= badflag
         endif

         if ( wetBulbTemperatureQua(j) .eq. 1 )  then           ! degC -> degK
               wetBulbTemperature(j)= wetBulbTemperature(j) +273.15
            else
               wetBulbTemperature(j)= badflag
         endif

         if ( precip24Hour(j) .eq. -99.9 )  then
               precip24Hour(j)= badflag
            else
               precip24Hour(j)= precip24Hour(j) *0.001   ! millimeter -> meter
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

c    -------     code figure transformed into visibility ( unit: m )  -------
      do j= 1,n
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

c  -----  code figure transformed into base of lowest cloud ( unit: m )  -----
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
         elseif ( skyLayerBase(1,j) .eq. 9. )  then
               skyLayerBase(1,j)= 3000.
         else
               skyLayerBase(1,j)= badflag
         endif

c ---- code figure transformed into base of cloud layer indicated (unit: m) ----
         if     ( skyLayerBase(2,j) .eq.  0. )  then
               skyLayerBase(2,j)= 15.
         elseif ( skyLayerBase(2,j) .gt.  0.  .and.
     ~            skyLayerBase(2,j) .lt. 51. )  then
               skyLayerBase(2,j)= skyLayerBase(2,j) *30.
         elseif ( skyLayerBase(2,j) .gt. 55.  .and.
     ~            skyLayerBase(2,j) .lt. 81. )  then
               skyLayerBase(2,j)= ( skyLayerBase(2,j) -50. ) *300.
         elseif ( skyLayerBase(2,j) .gt. 80.  .and.
     ~            skyLayerBase(2,j) .lt. 90. )  then
               skyLayerBase(2,j)= ( skyLayerBase(2,j) -74. ) *1500.
         elseif ( skyLayerBase(2,j) .eq. 90. )  then
               skyLayerBase(2,j)= 25.
         elseif ( skyLayerBase(2,j) .eq. 91. )  then
               skyLayerBase(2,j)= 75.
         elseif ( skyLayerBase(2,j) .eq. 92. )  then
               skyLayerBase(2,j)= 150.
         elseif ( skyLayerBase(2,j) .eq. 93. )  then
               skyLayerBase(2,j)= 250.
         elseif ( skyLayerBase(2,j) .eq. 94. )  then
               skyLayerBase(2,j)= 450.
         elseif ( skyLayerBase(2,j) .eq. 95. )  then
               skyLayerBase(2,j)= 800.
         elseif ( skyLayerBase(2,j) .eq. 96. )  then
               skyLayerBase(2,j)= 1250.
         elseif ( skyLayerBase(2,j) .eq. 97. )  then
               skyLayerBase(2,j)= 1750.
         elseif ( skyLayerBase(2,j) .eq. 98. )  then
               skyLayerBase(2,j)= 2250.
         elseif ( skyLayerBase(2,j) .eq. 99. )  then
               skyLayerBase(2,j)= 3000.
         else
               skyLayerBase(2,j)= badflag
         endif
      enddo

c --- code figure of total cloud cover or      transformed into metar format ---
c                    special cloud layer cover
      do 100 j= 1,n
      do 100 i= 1,maxSkyCover
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
100   continue

c               -------      dealing with lacking of data      -------
      do j= 1,n
         presWeather(j)= "UNK"

         dataPlatformType(j)= 1

         elevation(j)= 0.
         equivWindSpeed10m(j)= badflag
         precip1Hour(j)= badflag
         precip6Hour(j)= badflag 
         wmoId(j)= badflag 
      enddo

      go to 1000 

999   write (6,*) ' Error reading ship file'
      do j= 1,n
         write(6,*)  reportFlag(j),
     ~               stationName(j), latitude(j), longitude(j),
     ~               yy(j), mo(j), dd(j), hh(j), mn(j)
         write(6,*)  windDir(j), windSpeed(j), windQua(j),
     ~               visibility(j), presWeather(j),
     ~               seaLevelPress(j), seaLevelPressQua(j), 
     ~               temperature(j), temperatureQua(j),
     ~               skyCover(1,j), skyLayerBase(1,j)
         write(6,*)  dewpoint(j), dewpointQua(j), pressChangeChar(j),       
     ~               pressChange3Hour(j), pressChange3HourQua(j),
     ~               precip3Hour(j), wetBulbTemperature(j),
     ~               wetBulbTemperatureQua(j), windGust(j)
         write(6,*)  skyCover(2,j), skyLayerBase(2,j), precip24Hour(j)
      enddo

1000  return
      end
