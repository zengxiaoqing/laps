          subroutine write_aircraft_sub(lun,ext
     1                          ,a9_timeObs,a9_recptTime
     1                          ,i4time_sys
     1                          ,i4time_earliest          
     1                          ,i4time_latest            
     1                          ,latitude,longitude,altitude
     1                          ,windDir,windSpeed
     1                          ,temperature,relHumidity
     1                          ,l_geoalt
     1                          ,l_debug                           ! I
     1                          ,istat_ob)                         ! O

          character*(*) ext 

          character*9 a9_timeObs,a9_recptTime

          logical l_debug,l_geoalt

          real latitude, longitude

          istat_ob = 0

          call open_ext(lun,i4time_sys,ext(1:3),istatus)       

!         Test the altitude
          if(nanf(altitude) .eq. 1)then
              if(l_debug)write(6,*)' Altitude failed Nan test - reject'       
     1                            ,altitude
              goto 900
          endif

          if(altitude .gt. 20000. .or. altitude .lt. -1000.
     1                            .or. altitude .eq.     0.)then
              if(l_debug)write(6,*)' Altitude is suspect - reject'
     1                            ,altitude
              goto 900
          endif

!         Test the time
          call cv_asc_i4time(a9_timeObs,i4time_ob)
          if(i4time_ob .lt. i4time_earliest .OR.
     1       i4time_ob .gt. i4time_latest        )then ! outside time window
              if(l_debug)write(6,*)' time - reject '
     1           ,a9_timeObs,i4time_ob,i4time_earliest,i4time_latest
              goto 900        
          endif

          if(l_debug)write(6,1)a9_timeObs,a9_recptTime 
                     write(lun,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 
          istat_ob = 1

          if(l_geoalt)then
              if(l_debug)write(6,2)latitude,longitude,altitude
              write(lun,2)          latitude,longitude,altitude
 2            format(' Lat, lon, geoalt  '/f8.3,f10.3,f8.0)  
          else
              if(l_debug)write(6,3)latitude,longitude,altitude
              write(lun,3)          latitude,longitude,altitude
 3            format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  
          endif

!         Test for bad winds
!         if(char(dataDescriptor) .eq. 'X')then
!           if(char(errorType) .eq. 'W' .or. 
!    1         char(errorType) .eq. 'B'                         )then
!             if(l_debug)write(6,*)' QC flag is bad - reject wind'
!    1                 ,char(dataDescriptor),char(errorType)
!             goto 850
!           endif
!         endif

          if(abs(windSpeed) .gt. 250.)then
              if(l_debug)write(6,*)' wind speed is suspect - reject'
     1                              ,windSpeed

          elseif(int(windDir).lt.0 .or. int(windDir).gt.360)then     
              if(l_debug)write(6,*)' wind direction is suspect - reject'       
     1                              ,windDir

          else ! write out valid wind
              if(l_debug)then
                  write(6,4)int(windDir),windSpeed
              endif

              write(lun,4)int(windDir),windSpeed
 4            format(' Wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')     

          endif

!         Test/Write Temperature
          if(abs(temperature) .lt. 400.)then
              if(l_debug)write(6,13)temperature
                         write(lun,13)temperature
 13           format(' Temp:'/1x,f10.1)
       
          else
              if(l_debug)write(6,*)' Temperature is suspect - reject'
     1                             , temperature

          endif

!         Test/Write Relative Humidity
          if(relHumidity     .ge. 0.   .and. 
     1       relHumidity     .le. 1.00 
!    1       waterVaporQC(i) .le. 2    .and.
!    1       waterVaporQC(i) .ge. 0             
     1                                       )then
              if(l_debug)write(6,23)relHumidity
!             if(l_debug)write(6,*)' RH QC value = ',waterVaporQC(i)
              write(lun,23)relHumidity
 23           format(' RH:'/1x,f10.3)

          else
              if(l_debug)write(6,*)' RH rejected: '
     1                             ,relHumidity ! ,waterVaporQC(i)

          endif

 900      continue

!............................................................................

      return
      end
