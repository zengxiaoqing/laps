      subroutine read_radar_info(c4_radarname,rlat_radar,rlon_radar
     :                                       ,rheight_radar,istatus)

      character*4 c4_radarname
      character*80 c80_line

      call getenv('RADARNAME',c4_radarname)

      write(6,*)'read_radar_info: RADARNAME = ',c4_radarname

      open(11,file='radarinfo.dat',status='old')

 10   read(11,1,err=900)c80_line
 1    format(a)

      if(c80_line(1:4) .eq. c4_radarname)then ! We've got the right radar
          read(c80_line,2)ideg_lat,imin_lat,isec_lat
     1                   ,ideg_lon,imin_lon,isec_lon
     1                   ,iheight
 2        format(26x,i2,6x,i2,6x,i2,7x,i3,5x,i2,6x,i2,6x,i4)
          rlat_radar =  float(ideg_lat) + float(imin_lat)/60.
     1                                  + float(isec_lat)/3600.
          rlon_radar = -float(ideg_lon) - float(imin_lon)/60.
     1                                  - float(isec_lon)/3600.
          rheight_radar = iheight
          write(6,*)' rlat_radar,rlon_radar,rheight_radar '
     1               ,rlat_radar,rlon_radar,rheight_radar
          istatus = 1
          close(11)
          return
      endif       

      goto 10

 900  istatus = 0
      close(11)
      write(6,*)'read_radar_info: Radar data not found'
      return

 999  istatus = 1

      return
      end
