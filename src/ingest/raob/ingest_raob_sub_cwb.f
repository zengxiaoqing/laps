      subroutine get_raob_data_cwb ( i4time_sys, ilaps_cycle_time,
     ~             nx_l, ny_l, i4time_raob_earliest, i4time_raob_latest,
     ~             a9time_file, filename, istatus )

      integer   recNum  
      parameter ( recNum=80, levelNum=150 )

      character*(*)  filename
      character*9    a9time_file(recNum), a10_to_a9
      character*3    reportFlag(recNum)
      character*2    yy(recNum), mo(recNum), dd(recNum)
      character*2    hh(recNum), mn(recNum)
      character*10   time

      real  lat_a(nx_l,ny_l), lon_a(nx_l,ny_l), topo_a(nx_l,ny_l)
      real  elevation(recNum), latitude(recNum), longitude(recNum)
      real  pressure(recNum,levelNum), height(recNum,levelNum)
      real  temperature(recNum,levelNum)
      real  tempDewDiff(recNum,levelNum), dewpoint(recNum,levelNum)
      real  windDir(recNum,levelNum), windSpeed(recNum,levelNum)

      integer wmoId(recNum), logicRecNum, layerNum(recNum)
      integer heightQua(recNum,levelNum),temperatureQua(recNum,levelNum)
      integer dewpointQua(recNum,levelNum), windQua(recNum,levelNum)

      integer  d(12)
      data     d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif

      istatus= 0

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do i= 1,recNum
         read ( 1, 10, end=99, err=999 ) reportFlag(i), wmoId(i),
     ~               elevation(i), latitude(i), longitude(i),
     ~               yy(i), mo(i), dd(i), hh(i), mn(i), logicRecNum

         layerNum(i)= logicRecNum -2
         do 5 j= 1,layerNum(i)
5           read (1,20) pressure(i,j), height(i,j), heightQua(i,j),
     ~                  temperature(i,j), temperatureQua(i,j),        
     ~                  tempDewDiff(i,j), dewpointQua(i,j),
     ~                  windDir(i,j), windSpeed(i,j), windQua(i,j)

         read (1,*)

         if ( reportFlag(i) .ne. '*11' )  then
            write (6,*) 'read temp data heading error'
            go to 1000
         endif

         n= n+1
      enddo

10    format ( a3, i5, f4.0, 2f5.2, 2x, 5a2, i3 )
20    format ( 2x, f5.1, f5.0, i2, 2(f4.1,i2), 2f3.0, i2 )

c      ----------       examing data quality and changing units       ---------
99    do 100 i= 1,n      
         if ( elevation(i) .eq. -999. )  elevation(i)= r_missing_data

      do 100 j= 1,layerNum(i)

         if ( pressure(i,j) .eq. -999. )  pressure(i,j)= r_missing_data
         if ( heightQua(i,j) .ne. 1 )  height(i,j)= r_missing_data
         if ( temperatureQua(i,j).ne.1 ) temperature(i,j)=r_missing_data

         if ( temperatureQua(i,j).eq.1 .and. dewpointQua(i,j).eq.1 )then
               dewpoint(i,j)= temperature(i,j) -tempDewDiff(i,j)
            else
               dewpoint(i,j)= r_missing_data
         endif

         if ( windQua(i,j) .ne. 1 )  then
            windDir(i,j)= r_missing_data
            windSpeed(i,j)= r_missing_data
         endif

100   continue

c               ------ creat a9time_file in yydddhhmm format ------
      do i= 1,n
         if ( mn(i) .ne. '-9' )  then

            read (yy(i),'(i2)') iy
            read (mo(i),'(i2)') m1
            read (dd(i),'(i2)') id
            read (hh(i),'(i2)') ih
            read (mn(i),'(i2)') m2

            m2= m2 +30
            if ( m2 .ge. 60 )  then
               m2= m2 -60
               ih= ih +1

               if ( ih .ge. 24 )  then
                  ih= 0
                  id= id +1

                  if ( m1.eq.2  .and.  mod(iy,4).eq.0 )  d(m1)= d(m1) +1
                  if ( id .gt. d(m1) )  then
                     id= 1
                     m1= m1 +1
                         
                     if ( m1 .gt. 12 )  then
                        m1= 1
                        iy= iy +1
                     endif
                  endif
               endif
            endif
 
            call i2a ( iy, yy(i) )
            call i2a ( m1, mo(i) )
            call i2a ( id, dd(i) )
            call i2a ( ih, hh(i) )
            call i2a ( m2, mn(i) )

          else
            mn(i)= '00'
            if ( yy(i)(1:1) .eq. ' ' )  yy(i)= '0'//yy(i)(2:2)
            if ( mo(i)(1:1) .eq. ' ' )  mo(i)= '0'//mo(i)(2:2)
            if ( dd(i)(1:1) .eq. ' ' )  dd(i)= '0'//dd(i)(2:2)
            if ( hh(i)(1:1) .eq. ' ' )  hh(i)= '0'//hh(i)(2:2)
 
         endif
            
         time= yy(i)//mo(i)//dd(i)//hh(i)//mn(i)
         a9time_file(i)= a10_to_a9(time,istatus)
      enddo

      do 900 i= 1,n
         write (11,*) wmoId(i), elevation(i), latitude(i), longitude(i),
     ~                a9time_file(i), layerNum(i)

         do 900 j= 1,layerNum(i)
            write (11,*) height(i,j), pressure(i,j), temperature(i,j),
     ~                   dewpoint(i,j), windDir(i,j), windSpeed(i,j)
900   continue

      go to 1000

999   write (6,*) ' Error reading temp file'
      do i= 1,n
         write (6,*) reportFlag(i), wmoId(i),
     ~               elevation(i), latitude(i), longitude(i),
     ~               yy(i), mo(i), dd(i), hh(i), mn(i), layerNum(i)

         do j= 1,layerNum(i)
            write (6,*) pressure(i,j), height(i,j), heightQua(i,j),
     ~                  temperature(i,j), temperatureQua(i,j),
     ~                  tempDewDiff(i,j), dewpointQua(i,j),
     ~                  windDir(i,j), windSpeed(i,j), windQua(i,j)
         enddo
      enddo

1000  return
      end



      subroutine  i2a (ii,aa)

      character*2  aa
      integer      ii

      if ( ii .lt. 10 )  then
         write (aa,'(a1,i1)') '0',ii
      else
         write (aa,'(i2)') ii
      endif

      return
      end
