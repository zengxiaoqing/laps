
        subroutine get_domain_perimeter(ni,nj,LAPS_DOMAIN_FILE
     1                  ,lat,lon,topo
     1                  ,r_buffer,rnorth,south,east,west,istatus)

        real*4 lat(ni,nj),lon(ni,nj)
        real*4 topo(ni,nj)
        character*(*) LAPS_DOMAIN_FILE

        call get_laps_domain(ni,nj,LAPS_DOMAIN_FILE,lat,lon,topo
     1                      ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

        rnorth = -90.
        south  = +90.
        west = +1000.
        east = -1000.

        do i = 1,ni
        do j = 1,nj
            rnorth = max(rnorth,lat(i,j))
            south  = min(south ,lat(i,j))
            east   = max(east  ,lon(i,j))
            west   = min(west  ,lon(i,j))
        enddo ! j
        enddo ! i

        rnorth = rnorth + r_buffer
        south  = south  - r_buffer
        east   = east   + r_buffer
        west   = west   - r_buffer

        rnorth = min(rnorth, +90.)
        south  = max(south , -90.)
        east   = min(east  ,+180.)
        west   = max(west  ,-180.)

        write(6,101)rnorth,south,east,west
101     format(1x,' Lat/lon box around LAPS grid - NSEW ',4f9.2)

        return
        end

C -----------------------------------------------------------

        subroutine get_domain_perimeter_grid(ni,nj,LAPS_DOMAIN_FILE
     1                  ,lat,lon
     1                  ,r_buffer,rnorth,south,east,west,istatus)

        implicit none
        integer ni,nj,i,j
        integer istatus
        real east,west,rnorth,south
        real r_buffer
        real*4 lat(ni,nj),lon(ni,nj)
        character*(*) LAPS_DOMAIN_FILE
        logical lfnddlw,lfnddle
        logical lfndpln,lfndpls

c determine if either the e or w bndry cross the dateline
        lfnddle=.false.
        lfnddlw=.false.
        do j=2,nj
           if((lon(1,j-1).lt.0.0.and.lon(1,j).gt.0.0).or.
     &        (lon(1,j-1).gt.0.0.and.lon(1,j).lt.0.0))
     &lfnddlw=.true.
           if((lon(ni,j-1).lt.0.0.and.lon(ni,j).gt.0.0).or.
     &        (lon(ni,j-1).gt.0.0.and.lon(ni,j).lt.0.0))
     &lfnddle=.true. 
        enddo
c determine if either the n or s bndry cross the pole
        lfndpls=.false.
        lfndpln=.false.
        do i=2,ni
           if(lat(i-1,1).lt.0.0.and.lat(i,1).gt.0.0)
     &lfndpls=.true.
           if(lat(i-1,nj).gt.0.0.and.lat(i,nj).lt.0.0)
     &lfndpln=.true.
        enddo

c first pass: assume no dateline or pole boundaries
        rnorth = -90.
        south  = +90.
        west = +1000.
        east = -1000.
        do i = 1,ni
           rnorth = max(rnorth,lat(i,nj))
           south  = min(south ,lat(i,1))
        enddo
        do j = 1,nj
           east   = max(east  ,lon(ni,j))
           west   = min(west  ,lon(1,j))
        enddo ! j
c
c don't add buffer if too close to the dateline or pole boundary
        if(90.-rnorth.gt.r_buffer)then
           rnorth = rnorth + r_buffer
        endif
        if(south+90.0.gt.r_buffer)then
           south  = south  - r_buffer
        endif
        if(180.-east.gt.r_buffer)then
           east   = east   + r_buffer
        endif
        if(west+180.gt.r_buffer)then
           west   = west   - r_buffer
        endif

c do we have a boundary to deal with?
        
        if(lfnddlw)then      !dateline: want the min positive lon
           west = +1000.
           do j=1,nj
              if(lon(1,j).gt.0.0)
     &        west = min(west,lon(1,j))
           enddo
        endif

        if(lfnddle)then      !dateline: want the min negative lon
           east = +1000.
           do j=1,nj
              if(lon(ni,j).lt.0.0)
     &        east = min(east,lon(ni,j))
           enddo
        endif

        if(lfndpln)then
           rnorth=+90.
           do i=1,ni
              rnorth=min(rnorth,lat(i,nj))
           enddo
        endif

        if(lfndpls)then
           south=-90.
           do i=1,ni
              south=max(south,lat(i,1))
           enddo
        endif


c       rnorth = min(rnorth, +90.)
c       south  = max(south , -90.)
c       east   = min(east  ,+180.)
c       west   = max(west  ,-180.)

        return
        end

