      subroutine showgrid(lat,lon,nx,ny,delta,c6_proj,std_lat,std_lat2,
     +     std_lon, cenlat,cenlon)
      implicit none
      integer nx,ny, i, j
      real lat(nx,ny), lon(nx,ny)
      real ll_lat, ll_lon, ur_lat,ur_lon
      real delta,std_lat,std_lat2,std_lon,cenlat,cenlon
      character*6 c6_proj

      integer proj_idx
      integer cornertype
      integer istatus

      if(c6_proj.eq.'plrstr') then
         proj_idx=1
         std_lat2=0.
      else if(c6_proj.eq.'lambrt') then
         proj_idx=3
      else if(c6_proj.eq.'merctr') then
         proj_idx=9
         std_lat2=0.
      else
         print*,c6_proj,' is not a recognized projection'
         return
      endif
      
c
c we want a background slightly larger than the laps grid
c      
      ll_lat = lat(1,1)
      ll_lon = lon(1,1)
      ur_lat = lat(nx,ny)
      ur_lon = lon(nx,ny)

c      print *,ll_lat,ll_lon,ur_lat,ur_lon
c      print *,lon(1,1),lon(nx,1),lon(1,ny),lon(nx,ny)
      do j=1,ny
         ll_lon=min(ll_lon,lon(1,j))
         ur_lon=max(ur_lon,lon(nx,j))
      enddo
      do i=1,nx
         ll_lat=min(ll_lat,lat(i,1))
         ur_lat=max(ur_lat,lat(i,ny))
      enddo

c      print *,ll_lat,ll_lon,ur_lat,ur_lon

      ll_lat = ll_lat - (lat(3,3)-lat(1,1))
      ll_lon = ll_lon - (lon(3,3)-lon(1,1))


      ur_lat = ur_lat + (lat(3,3)-lat(1,1))
      ur_lon = ur_lon + (lon(3,3)-lon(1,1))


c      print *,ll_lat,ll_lon,ur_lat,ur_lon

      call create_bcd_bkgnd(proj_idx,std_lat,std_lon,std_lat2,ll_lat
     +     ,ll_lon,ur_lat,ur_lon,2,istatus)

      open(2,file="gridpoints")
      open(3,file="centerpoint")
      write(3,'(f16.8,a,f16.8)') cenlon,',',cenlat
      close(3)
c      do j=1,ny
c         do i=1,nx
c            write(2,'(2f16.8)') lon(i,j),lat(i,j)
c         enddo
c      enddo

      if(lon(1,1).eq.lon(1,ny)) then
         write(2,103) lon(1,1),lat(1,1),lon(nx,ny),lat(nx,ny)
      else
         do j=1,ny-1,5
            write(2,103) lon(1,j),lat(1,j)
     +           ,lon(1,MIN(j+5,NY)),lat(1,MIN(j+5,NY))
         enddo
         
      endif
      if(lon(nx,1).eq.lon(nx,ny)) then
         write(2,103) lon(nx,1),lat(nx,1),lon(nx,ny),lat(nx,ny)
      else
         do j=1,ny-1,5
            write(2,103) lon(nx,j),lat(nx,j)
     +           ,lon(nx,MIN(j+5,NY)),lat(nx,MIN(j+5,NY))
         enddo
         
      endif

      if(lat(1,1).eq.lat(nx,1)) then
         write(2,103) lon(1,1),lat(1,1),lon(nx,1),lat(nx,1)
      else
         do i=1,nx-1,5
            write(2,103) lon(i,1),lat(i,1)
     +           ,lon(MIN(i+5,NX),1),lat(MIN(i+5,NX),1)
         enddo
         
      endif
      if(lat(1,ny).eq.lat(nx,ny)) then
         write(2,103) lon(1,ny),lat(1,ny),lon(nx,ny),lat(nx,ny)
      else
         do i=1,nx-1,5
            write(2,103) lon(i,ny),lat(i,ny)
     +           ,lon(MIN(i+5,NX),ny),lat(MIN(i+5,NX),ny)
         enddo
         
      endif

c      do i=2,nx
c         write(2,103) lon(i-1,1),lat(i-1,1),lon(i,1),lat(i,1)
c         write(2,103) lon(i-1,ny),lat(i-1,ny),lon(i,ny),lat(i,ny)
c      enddo

      write(2,104) ll_lon,ur_lon
      write(2,105) ll_lat,ur_lat
c      write(2,*) 'plot "centerpoint" ls 3' 

 103  format('set arrow from ',f16.8,',',f16.8,' to '
     +     ,f16.8,',',f16.8,' nohead')

 104  format('set xrange [',f16.8,':',f16.8,']')
 105  format('set yrange [',f16.8,':',f16.8,']')

      close(2)

      return
      end

      subroutine get_gridnl(mode)
      implicit none
      integer mode
      integer len
      integer istatus
      character*256 directory
      character*256 fname
      character*200 cdataroot
      character*10  c10_grid_fname
      namelist /grid_nl/ mode

      mode = 0
      call find_domain_name(cdataroot,c10_grid_fname,istatus)
      if(istatus.ne.1)then
         print*,'Error returned from find_domain_name'
         return
      endif
      call get_directory(c10_grid_fname,directory,len)
      fname = directory(1:len)//'grid.nl'
      open(3,file=fname,status='old',err=101)
      read(3,grid_nl,err=101,end=101)
      close(3)

 101  return
      end
