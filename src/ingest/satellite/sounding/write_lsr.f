      subroutine write_lsr(nxl,nyl,nch,csatid,wavelength,
     &data,i4time,istatus)

      implicit none

      integer nxl,nyl
      integer nch

      integer i4time,i4time_nearest
      integer istatus
      integer i,j,k
      integer ni
      integer len_lsr
      integer lvl_lsr(nch)
      integer izero

      real   data(nxl,nyl,nch)
      real   wavelength(nch)

      character*125 c_lsr(nch)
      character*10 units_lsr(nch)
      character*3 var_lsr(nch)
      character*4 lvl_coord_lsr(nch)
      character*150 dir_lsr
      character*31 ext_lsr
      character*6  csatid
      character*2  c_num
      character*4  cw
c ---------------------------------------
c Output for LAPS lsr files as indicated.
c
      call get_directory('lsr',dir_lsr,len_lsr)
      ni=index(csatid,' ')-1
      if(ni.le.0)ni=6
      dir_lsr=dir_lsr(1:len_lsr)//csatid(1:ni)//'/'
      ext_lsr = 'lsr'
c
      do i=1,nch
         lvl_lsr(i) = 0
         lvl_coord_lsr(i) = 'AGL'
         units_lsr(i) = 'Radiance'
         c_num = ' '
         write(c_num,100)i 
100      format(i2)

         if(c_num(1:1).eq.' ')then
            c_num(1:1)='0'
         endif
         var_lsr(i) = 's'//c_num

         write(cw,200)wavelength(i)
200      format(f4.1)
         do j=1,4
            if(cw(j:j).eq.' ')cw(j:j)='0'
         enddo

         c_lsr(i)=csatid//' SAT sounding radiance: '//cw//' (u)'

      enddo !i

      print*,'checking for <= zero radiance: before write'
      do k=1,nch
         call zero_rad_check(nxl,nyl,k,data(1,1,k))
      enddo

      call write_laps_data(i4time,
     &                     dir_lsr,
     &                     ext_lsr,
     &                     nxl,nyl,
     &                     nch,
     &                     nch,
     &                     var_lsr,
     &                     lvl_lsr,
     &                     lvl_coord_lsr,
     &                     units_lsr,
     &                     c_lsr,
     &                     data,
     &                     istatus)

      if(istatus.ne.1)then
         write(6,*)'Error writing lsr - write_laps_data'
      endif


      if(.false.)then

       print*,'spot check for <= zero radiance'

       do k=1,nch

         call get_2dgrid_dname(dir_lsr
     1         ,i4time,900,i4time_nearest
     1         ,ext_lsr,var_lsr(k),units_lsr(k)
     1         ,c_lsr,nxl,nyl,data(1,1,k),0,istatus)

         if(istatus.ne.1)then
            print*,'subroutine error get_2dgrid_dname'
            return
         endif
   

         print*,'lsr 2d field obtained.'
         call zero_rad_check(nxl,nyl,k,data(1,1,k))
         print*,'finished spot check'

       enddo

      endif

      return
      end

c
      subroutine zero_rad_check(nxl,nyl,k,data)

      real   data(nxl,nyl)

      print*,'checking for <= zero value radiance'
      print*,'Channel ',k
      izero=0
      do j=1,nyl
      do i=1,nxl
         if(data(i,j).le.0.0)then
            print*,'i/j/2d: ',i,j,data(i,j)
            izero=izero+1
         endif
      enddo
      enddo
      if(izero.gt.0)print*,'found ',izero,' <= zero radiance'

      return
      end
