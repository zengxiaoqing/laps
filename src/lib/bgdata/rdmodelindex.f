      subroutine readindexfile(cfname,nfldsmax,nlevsmax
     +,nflds,ilevs,p,ivarcoord,ivarid,istatus)
c
c 7-98: J. Smart.   routine reads yyjjjhhmmhhmm.index files for AFWA
c                   models (AVN, NOGAPS, etc). Routine is used in readdgprep.f
c                   and readnogaps. Result is used to find appropriate model
c                   fields in AFWA degrib files.
c
      implicit none

      integer  nfldsmax
      integer  nlevsmax
      integer  nflds
      integer  ivid
      integer  ivrcd
      integer  ilevs(nfldsmax)
      integer  ivarid(nfldsmax)
      integer  ivarcoord(nfldsmax)

      integer  i,j,ip,ifl
      integer  lun
      integer  istatus
      integer  iostatus

      real     p(nlevsmax,nfldsmax)

      character cfname*(*)

      istatus=1  !default error return (1=failure)
      lun=120
      ifl=index(cfname,' ')-1

      open(lun,file=cfname,form='formatted',iostat=iostatus,
     +     status='old',err=500)


      nflds=0
      do i=1,nfldsmax
         read(lun,*,end=502,err=501)ivid
         nflds=nflds+1
         ivarid(nflds)=ivid
         read(lun,*,end=502,err=501)ivrcd
         ivarcoord(nflds)=ivrcd
         read(lun,*,end=502,err=501)ilevs(nflds)
         do j=1,ilevs(i)
            read(lun,*,end=502,err=501)ip
            p(j,nflds)=float(ip)
         enddo
      enddo

502   close(lun)
      istatus=0  !success reading index file
      return

500   print*,'error opening file ',cfname(1:ifl)
      print*,'iostatus = ',iostatus
      return
501   print*,'error reading file ',cfname(1:ifl)
      return

C502   print*,'end of file encountered'

      end
