      subroutine read_buoy_cwb ( filename, recNum, 
     +     dataPlatformType, dewpoint, elevation, equivWindSpeed10m,
     +     latitude, longitude, precip1Hour, precip24Hour,
     +     precip6Hour, presWeather, pressChange3Hour,
     +     pressChangeChar, seaLevelPress, seaSurfaceTemp,
     +     stationName, temperature, timeObs,
     +     visibility, wetBulbTemperature, windDir, windGust,
     +     windSpeed, badflag, istatus)
 
      integer recNum, i4time

      character*25 presWeather(recNum)
      character*8  stationName(recNum)
      integer  dataPlatformType(recNum), pressChangeChar(recNum)

      real  dewpoint(recNum), elevation(recNum),
     ~      equivWindSpeed10m(recNum), latitude(recNum),
     ~      longitude(recNum), precip1Hour(recNum),
     ~      precip24Hour(recNum), precip6Hour(recNum),
     ~      pressChange3Hour(recNum), seaLevelPress(recNum),
     ~      seaSurfaceTemp(recNum), temperature(recNum),
     ~      visibility(recNum), wetBulbTemperature(recNum),
     ~      windDir(recNum), windGust(recNum), windSpeed(recNum)

      double precision  timeObs(recNum)

      character*(*) filename
      character*3   reportFlag(recNum)
      character*2   yy(recNum), mo(recNum), dd(recNum)
      character*2   hh(recNum), mn(recNum)
      character*10  time(recNum)
      character*9   a10_to_a9

      integer  logicRecNum(recNum)
      integer  windQua(recNum), temperatureQua(recNum)
      integer  seaLevelPressQua(recNum), seaSurfaceTempQua(recNum)
      integer  dewpointQua(recNum)

      istatus = 0
 
      open ( 1, file=filename, status='old' )

      n= 0
      do j= 1,recNum
         read ( 1, 10, end=99, err=999 ) reportFlag(j),
     ~               stationName(j), latitude(j), longitude(j),
     ~               yy(j), mo(j), dd(j), hh(j), mn(j), logicRecNum(j)
         read (1,20) windDir(j), windSpeed(j), windQua(j),
     ~               temperature(j), temperatureQua(j),
     ~               seaLevelPress(j), seaLevelPressQua(j),
     ~               pressChange3Hour(j)
         read (1,30) seaSurfaceTemp(j), seaSurfaceTempQua(j),
     ~               dewpoint(j), dewpointQua(j)

         do 5 i= 1,logicRecNum(j)-3 
5           read (1,*)

         if ( reportFlag(j) .ne. '*81' )  then
            write (6,*) 'read data heading error'
            go to 1000
         endif

         n= n+1
      enddo

10    format ( a3, a5, 4x, 2f5.2, 2x, 5a2, i3 )
20    format ( 2f3.0, i1, f4.1, i1, f5.1, i1, x, f3.1 ) 
30    format ( f4.1, i1, 22x, f5.1, i1 )

c      ----------       examing data quality and changing units       ---------
99    do j= 1,n
         if ( windQua(j) .eq. 9 )  then
            windDir(j)= badflag
            windSpeed(j)= badflag
         endif

         if ( temperatureQua(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15         ! degC -> degK
            else
               temperature(j)= badflag
         endif

         if ( seaLevelPressQua(j) .eq. 1 )  then
               seaLevelPress(j)= seaLevelPress(j) *100        ! mb -> pascal
            else
               seaLevelPress(j)= badflag
         endif
  
         if ( seaSurfaceTempQua(j) .eq. 1 )  then
               seaSurfaceTemp(j)= seaSurfaceTemp(j) +273.15   ! degC -> degK
            else
               seaSurfaceTemp(j)= badflag
         endif

         if ( dewpointQua(j) .eq. 1 )  then
               dewpoint(j)= dewpoint(j) +273.15               ! degC -> degK
            else
               dewpoint(j)= badflag
         endif

         if ( pressChange3Hour(j).eq.-9.9 ) pressChange3Hour(j)= badflag        

         elevation(j)= 0

         if ( yy(j)(1:1) .eq. ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( mo(j)(1:1) .eq. ' ' )  mo(j)= '0'//mo(j)(2:2)
         if ( dd(j)(1:1) .eq. ' ' )  dd(j)= '0'//dd(j)(2:2)
         if ( hh(j)(1:1) .eq. ' ' )  hh(j)= '0'//hh(j)(2:2)
         if ( mn(j)(1:1) .eq. ' ' )  mn(j)= '0'//mn(j)(2:2)
         time(j)= yy(j)//mo(j)//dd(j)//hh(j)//mn(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeObs(j)= dble( i4time )                    ! seconds since 1960
      enddo

c               -------      dealing with lacking of data      -------
      do j= 1,n
         presWeather(j)= "UNK"

         pressChangeChar(j)= int(badflag)
         dataPlatformType(j)= int(badflag)

         equivWindSpeed10m(j)= badflag
         precip1Hour(j)= badflag
         precip24Hour(j)= badflag
         precip6Hour(j)= badflag 
         pressChangeChar(j)= badflag
         visibility(j)= badflag
         wetBulbTemperature(j)= badflag
         windGust(j)= badflag
      enddo

      istatus= 1
*     go to 1000 

999   do j= 1,n
         write(6,*) reportFlag(j),
     ~              stationName(j),latitude(j), longitude(j),
     ~              yy(j), mo(j), dd(j), hh(j), mn(j), logicRecNum(j)
         write(6,*) windDir(j), windSpeed(j), windQua(j),
     ~              temperature(j), temperatureQua(j),
     ~              seaLevelPress(j), seaLevelPressQua(j),
     ~              pressChange3Hour(j)
         write(6,*) seaSurfaceTemp(j), seaSurfaceTempQua(j),
     ~              dewpoint(j), dewpointQua(j), timeobs(j)
      enddo

1000  return
      end
