      subroutine write_lsr(nxl,nyl,nch,csatid,wavelength,
     &data,i4time,istatus)

      implicit none

      integer*4 nxl,nyl
      integer*4 nch

      integer*4 i4time
      integer*4 istatus
      integer*4 i,j
      integer*4 len_lsr

      real*4 data(nxl,nyl,nch)
      real*4 wavelength(nch)

      character*125 c_lsr(nch)
      character*10 units_lsr(nch)
      character*3 var_lsr(nch)
      character*4 lvl_coord_lsr(nch)
      integer*4 lvl_lsr(nch)
      character*150 dir_lsr
      character*31 ext_lsr
      character*5  csatid
      character*2  c_num
      character*4  cw
c ---------------------------------------
c Output for LAPS lsr files as indicated.
c
      call get_directory('lsr',dir_lsr,len_lsr)
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

      return
      end
