      subroutine readvrclut(cradtype,nx,ny,ri,rj,istatus)
c
c

      implicit none

ccc      include 'lapsparms.for'

      integer*4 nx,ny
      integer*4 istatus
      integer*4 n
      integer*4 i,n1,n2

      real*4    ri(nx,ny)
      real*4    rj(nx,ny)
      real*4    rdummy(nx,ny)

      character*20  cnm
      character*50  path
      character*255 table_path
      character*3   cradtype

      istatus = 0
      call get_directory('static',path,n1)
      path=path(1:n1)//'vrc/ '
      n1=index(path,' ')-1
      cnm='wsi_llij_lut_'
      n2=index(cnm,' ')-1

ccc      do i = 1,n_radar_types

ccc        if(cradtype.eq.c_raddat_types(i))then

         table_path=path(1:n1)//cnm(1:n2)//cradtype//'.lut'
         n=index(table_path,' ')

         call read_table(table_path,nx,ny,
     &        rdummy,rdummy,ri,rj,istatus)

ccc         endif
ccc      enddo

      goto 1000

23    write(6,*)'Error reading or eof ll/ij lookup table'
      istatus = -1
      goto 1000

101   write(6,*)'Error openning file ',table_path(1:n)
      istatus = -1

1000  return
      end
