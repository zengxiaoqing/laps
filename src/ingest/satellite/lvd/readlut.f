      subroutine readlut(csat_id,csat_typ,maxch,nft,
     &ntm,chtype,nx,ny,ri,rj,istatus)
c
      implicit none

      integer nx,ny
      integer maxch
      integer istatus
      integer i,j,n,np
      integer lend,lenf
      integer ispec
      integer istat
      integer nft
      integer ntm(nft)

      logical   lgot_lut(maxch)

      real*4    ri(nx,ny,maxch)
      real*4    rj(nx,ny,maxch)
      real*4    ri_in(nx,ny)
      real*4    rj_in(nx,ny)
      real*4    rdummy(nx,ny)

      character*100 cpath
      character*255 file
      character*6   csat_id
      character*3   csat_typ
      character*3   chtype(maxch,nft)
      character*3   ct
c
c-----------------------------------------------------------------------
c
      istatus = -1

      call get_directory('static',cpath,lend)
      cpath=cpath(1:lend)//'lvd/'
      lend=index(cpath,' ')-1

      do i=1,maxch
         lgot_lut(i)=.false.
      enddo

      do j=1,nft
      do i=1,ntm(j)

         call lvd_file_specifier(chtype(i,j),ispec,istat)
         if(istat.eq.0)then

            if(.not.lgot_lut(ispec))then
              lgot_lut(i)=.true.

              ct=chtype(i,j)
              goto(3,2,3,2,2)ispec
2             ct='ir'

3             n=index(ct,' ')-1
              if(n.le.0)n=3
              file=cpath(1:lend)//csat_id//'-llij-'
              lenf=index(file,' ')-1
              file=file(1:lenf)//ct(1:n)//'-'//csat_typ//'.lut'

              n=index(file,' ')-1
              open(12,file=file,
     &form='unformatted',status='old',err=101)
              write(6,*)'Reading ',file(1:n)
              read(12,err=23,end=23) rdummy
              read(12,err=23,end=23) rdummy
              read(12,err=23,end=23) ri_in
              read(12,err=23,end=23) rj_in
              close (12)

              call move(ri_in,ri(1,1,ispec),nx,ny)
              call move(rj_in,rj(1,1,ispec),nx,ny)
           endif
         endif
      enddo
      enddo

      istatus = 0
      goto 1000

23    write(6,*)'Error reading or eof ll/ij lookup table'
      goto 1000

101   write(6,*)'Error opening file ',file(1:n)

1000  return
      end
