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
      subroutine read_buoy_cwb ( filename ,maxSkyCover, recNum,
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

      integer  logicRecNum(recNum)
      integer  windQua(recNum), temperatureQua(recNum)
      integer  seaLevelPressQua(recNum), seaSurfaceTempQua(recNum)
      integer  dewpointQua(recNum)

      real  relaHumility(recNum)

      istatus = 0
      n= 0
 
      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do j= 1,recNum
         read ( 1, 10, end=99, err=999 )
     ~               reportFlag(j), wmoId(j), latitude(j), longitude(j),       
     ~               yy(j), mo(j), dd(j), hh(j), mn(j), logicRecNum(j)
         read (1,20) windDir(j), windSpeed(j), windQua(j),
     ~               temperature(j), temperatureQua(j),
     ~               seaLevelPress(j), seaLevelPressQua(j),
     ~               pressChangeChar(j), pressChange3Hour(j)
         read (1,30) seaSurfaceTemp(j), seaSurfaceTempQua(j),
     ~               dewpoint(j), dewpointQua(j)
         read (1,40) relaHumility(j)

         do 5 i= 1,logicRecNum(j)-4
5           read (1,*)

         if ( reportFlag(j) .ne. '*81' )  then
            write (6,*) 'read buoy data heading error'
            go to 1000
         endif

         n= n+1
      enddo

10    format ( a3, i5, 4x, 2f5.2, 2x, 5a2, i3 )
!HJ: W>=D+3. 10/14/2013
20    format ( 2f3.0, i1, f4.1, i1, f5.1, 2i1, f4.1 ) 
30    format ( f4.1, i1, 22x, f5.1, i1 )
40    format ( 30x, f3.0 )

c      ----------       examing data quality and changing units       ---------
99    do j= 1,n
         if ( windQua(j) .eq. 9 )  then
            windDir(j)= badflag
            windSpeed(j)= badflag
         endif

         if ( seaLevelPressQua(j) .eq. 1 )  then
               seaLevelPress(j)= seaLevelPress(j) *100.         ! mb -> pascal
            else
               seaLevelPress(j)= badflag
         endif
  
         if ( pressChange3Hour(j) .eq. 1 )  then
               pressChange3Hour(j)= pressChange3Hour(j) *100.   ! mb -> pascal
            else
               pressChange3Hour(j)= badflag
         endif

         if ( seaSurfaceTempQua(j) .eq. 1 )  then
               seaSurfaceTemp(j)= seaSurfaceTemp(j) +273.15     ! degC -> degK
            else
               seaSurfaceTemp(j)= badflag
         endif

         if ( dewpointQua(j) .eq. 1 )  then
               dewpoint(j)= dewpoint(j) +273.15                 ! degC -> degK
            else
               if ( temperature(j) .ne. -999.9   .and. 
     ~              relaHumility(j) .ne. -99. )  then
                  dewpoint(j)= dwpt( temperature(j), relaHumility(j) )
                  dewpoint(j)= dewpoint(j) +273.15              ! degC -> degK
                else
                  dewpoint(j)= badflag
                  relaHumility(j)= badflag
               endif
         endif

         if ( temperatureQua(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15           ! degC -> degK
            else
               temperature(j)= badflag
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

c               -------      dealing with lacking of data      -------
      do j= 1,n
         presWeather(j)= "UNK"
         skyCover(1,j)=  "UNK"
         skyCover(2,j)=  "UNK"
         stationName(j)= "UNK"

         pressChangeChar(j)= int(badflag)
         dataPlatformType(j)= 0

         elevation(j)= 0.
         equivWindSpeed10m(j)= badflag
         precip1Hour(j)= badflag
         precip24Hour(j)= badflag
         precip3Hour(j)= badflag 
         precip6Hour(j)= badflag 
         skyLayerBase(1,j)= badflag 
         skyLayerBase(2,j)= badflag 
         visibility(j)= badflag
         wetBulbTemperature(j)= badflag
         windGust(j)= badflag
      enddo

      go to 1000 

999   write (6,*) ' Error reading buoy file'
      do j= 1,n
         write(6,*) reportFlag(j), wmoId(j),latitude(j), longitude(j),
     ~              yy(j), mo(j), dd(j), hh(j), mn(j), logicRecNum(j)
         write(6,*) windDir(j), windSpeed(j), windQua(j),
     ~              temperature(j), temperatureQua(j),
     ~              seaLevelPress(j), seaLevelPressQua(j),
     ~              pressChangeChar(j), pressChange3Hour(j)
         write(6,*) seaSurfaceTemp(j), seaSurfaceTempQua(j),
     ~              dewpoint(j), dewpointQua(j), timeobs(j)
         write(6,*) relaHumility(j)
      enddo

1000  return
      end
