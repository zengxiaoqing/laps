      subroutine get_cloud_drift_cwb_hdsw 
     ~           (i4time_sys, i4_window, nx_l, ny_l, filename, istatus)

      parameter  ( loopNum=500 )

      character*(*)  filename
      character*2    yy, mo, dd, hh, mn
      character*9    a9timeObs(loopNum), a10_to_a9
      character*10   time

      real  lat_a(nx_l,ny_l), lon_a(nx_l,ny_l), topo_a(nx_l,ny_l)
      real  latitude(loopNum), longitude(loopNum)
      real  pressure(loopNum), windDir(loopNum), windSpeed(loopNum)

      integer  pressureQua(loopNum), windQua(loopNum), recNum, inNum

      call get_r_missing_data (r_missing_data,istatus)
      if ( istatus .ne. 1 )  then
         write (6,*) 'Error getting r_missing_data'
         return
      endif

      recNum= 0
      inNum= 0              !  inNum : the record number within time window
      istatus= 0

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do i= 1,loopNum
         read (1,50,end=99,err=9) latitude(i), longitude(i), 
     ~              yy, mo, dd, hh, mn, pressure(i), pressureQua(i),
     ~              windDir(i), windSpeed(i), windQua(i)

c               ------ creat a9timeObs in yydddhhmm format ------
         if ( yy(1:1) .eq. ' ' )  yy= '0'//yy(2:2)
         if ( mo(1:1) .eq. ' ' )  mo= '0'//mo(2:2)
         if ( dd(1:1) .eq. ' ' )  dd= '0'//dd(2:2)
         if ( hh(1:1) .eq. ' ' )  hh= '0'//hh(2:2)
         if ( mn(1:1) .eq. ' ' )  mn= '0'//mn(2:2)
 
         time= yy//mo//dd//hh//mn
         a9timeObs(i)= a10_to_a9(time,istatus)
	 if ( istatus .ne. 1 )  then
	    write (6,*) 'Bad observation time - reject ', time
	 endif  

         call cv_asc_i4time( a9timeObs(i), i4time_obs )
         i4_residue= abs( i4time_obs -i4time_sys )

c          ----------    test if raob is within time window    ----------
	 if ( i4_residue .le. i4_window )  then
	    inNum= inNum +1

	    latitude(inNum)= latitude(i)
	    longitude(inNum)= longitude(i)
	    pressure(inNum)= pressure(i)
	    pressureQua(inNum)= pressureQua(i)
            windDir(inNum)= windDir(i)
	    windSpeed(inNum)= windSpeed(i)
	    windQua(inNum)= windQua(i)
	    a9timeObs(inNum)= a9timeObs(i)
         endif
         go to 10

9        write (6,*) 'Error reading variables of hdsw code '
         write (6,*) latitude(i), longitude(i),
     ~               yy, mo, dd, hh, mn, pressure(i), pressureQua(i),
     ~               windDir(i), windSpeed(i), windQua(i)

10       recNum= recNum +1       
      enddo

50    format ( 11x, 2f5.2, 5a2, f4.0, i1, f3.0, f4.1, i1 )

c      ----------       examine data quality and change units       ---------
99    do 100 i= 1,inNum
         if ( pressureQua(i) .eq. 1 )  then
	    pressure(i)= pressure(i) *100.          !  unit: mb --> pa
          else
	    pressure(i)= r_missing_data
         endif

         if ( windQua(i) .ne. 1 )  then
            windDir(i)= r_missing_data
            windSpeed(i)= r_missing_data
         endif
100   enddo

      do 900 i= 1,inNum
         call open_ext(11,i4time_sys,'cdw',istatus)
900      write (11,21) latitude(i), longitude(i), pressure(i), 
     ~                 windDir(i), windSpeed(i), a9timeObs(i)
 21      format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

      write (6,*) 'found', inNum,'high density satellite wind data ',
     ~            'within time window in ', filename
       
1000  return
      end
