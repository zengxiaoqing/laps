      subroutine read_gvap (filename, nstations, lat,lon, wt,w1,w2,w3,
     1 nn, istatus)

c     routine to read goes watervapor for system testing
c     author Dan Birkenheuer
c     2/8/99
c     
c     paramter list
c     filename, name of the wisconsin file
c     nstations, dimension of variables
c     lat, latitude of each station
c     lon, longitude of each station
c     wt, total water
c     w1, layer 1 water
c     w2, layer 2 water
c     w3, layer 3 water
c     nn, number of stations with real data for consideration
c
      implicit none

c     input variables
      character*9 filename
      integer nstations,nn,istatus,idummy
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)

c     internal variables

      integer i
      real dummy

c     reading goes 8
      
      open(22, file='/data/rapb/taiwan/goeswv/'//filename//'.wv8',form=
     1 'formatted',status='old',err = 668)
      read(22,*,end=668,err=668)   ! first header line is ignored
      do i = 1,nstations
         read(22,*,end=665,err=665) idummy,lat(i),lon(i),
     1              idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo

 665  close (22)
      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read GOES 8'
      istatus = 1
c      write(6,*) (wt(i),i=1,nn)

      go to 669

 668  write(6,*) 'failed reading GOES 8'
      nn = 1

 669  continue

c     reading goes 10

      filename = filename(1:7)//'24'

      open (23, file='/data/rapb/taiwan/goeswv/'//filename//'.wv10',
     1     form='formatted',status='old',err = 666)
      read(23,*,end=666,err=666)   ! first header line is ignored
      do i = nn,nstations
         read(23,*,end=667,err=667) idummy,lat(i),lon(i),
     1              idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo


 667  close(23)

      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read GOES 10'
      istatus = 1
c      write(6,*) (wt(i),i=1,nn)

      return

 666  if(nn.eq.0) then
         istatus = 0
         write(6,*) ' no available gvap data'
      else
         write(6,*) 'GOES 10 not avail but GOES 8 data detetected'
         istatus = 1
      endif

      return

      end

