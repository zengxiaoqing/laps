
        subroutine get_latlon_perimeter(ni,nj,r_buffer             ! I
     1                  ,lat,lon,topo                              ! O
     1                  ,rnorth,south,east,west,istatus)           ! O

cdoc    Obtain lat/lon box surrounding laps grid, including a buffer
cdoc    Works on the basis of max/min lat/lons and limited accounting
cdoc    for when the dateline or poles are near the box.
cdoc    This version is more generic with static file source.

        real lat(ni,nj),lon(ni,nj)
        real topo(ni,nj),fracland(ni,nj)

        call get_laps_domain_95(ni,nj,lat,lon,topo,fracland,gridsp
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
101     format(1x,' Lat/lon box around LAPS grid - NSEW ',4f10.3)

        return
        end

        subroutine get_domain_perimeter(ni,nj,LAPS_DOMAIN_FILE
     1                  ,lat,lon,topo
     1                  ,r_buffer,rnorth,south,east,west,istatus)

cdoc    Obtain lat/lon box surrounding laps grid, including a buffer
cdoc    Works on the basis of max/min lat/lons and limited accounting
cdoc    for when the dateline or poles are near the box.

        real lat(ni,nj),lon(ni,nj)
        real topo(ni,nj)
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
101     format(1x,' Lat/lon box around LAPS grid - NSEW ',4f10.3)

        return
        end

C -----------------------------------------------------------

        subroutine get_domain_perimeter_grid(ni,nj,LAPS_DOMAIN_FILE
     1                  ,lat,lon
     1                  ,r_buffer,rnorth,south,east,west,istatus)

cdoc    Obtain lat/lon box surrounding laps grid, including a buffer
cdoc    Works on the basis of lat/lons and greater accounting
cdoc    for when the dateline or poles are near or in the box.

        implicit none
        integer ni,nj,i,j
        integer istatus
        real east,west,rnorth,south
        real r_buffer
        real lat(ni,nj),lon(ni,nj)
        character*(*) LAPS_DOMAIN_FILE
        logical lfnddateline
        logical lfndgrenwich
        logical lfndpolen
        logical lfndpoles

c determine if either the e or w bndry cross the dateline
        lfnddateline=.false.
        lfndgrenwich=.false.
        do j=2,nj
           if((lon(1,j-1).lt.0.0.and.lon(1,j).ge.0.0).or.
     &        (lon(1,j-1).ge.0.0.and.lon(1,j).lt.0.0))
     &         lfnddateline=.true.
           if((lon(ni,j-1).lt.0.0.and.lon(ni,j).ge.0.0).or.
     &        (lon(ni,j-1).ge.0.0.and.lon(ni,j).lt.0.0))
     &         lfndgrenwich=.true. 
        enddo

        do i=2,ni
           if((lon(i-1,1).lt.0.0.and.lon(i,1).gt.0.0).or.
     &        (lon(i-1,nj).lt.0.0.and.lon(i,nj).gt.0.0))
     &         lfndgrenwich=.true.
           if((lon(i-1,1).gt.0.0.and.lon(i,1).lt.0.0).or.
     &        (lon(i-1,nj).gt.0.0.and.lon(i,nj).lt.0.0))
     &         lfnddateline=.true. 
        enddo

        lfndpoles=.false.
        lfndpolen=.false.
c cross the north pole?
        do j=2,nj
           if((lat(1,j-1).gt.0.0.and.lat(1,j).lt.0.0).or.
     &        (lat(ni,j-1).gt.0.0.and.lat(ni,j).lt.0.0))
     &        lfndpolen=.true.
c cross the equator?
           if((lat(1,j-1).lt.0.0.and.lat(1,j).ge.0.0).or.
     &        (lat(ni,j-1).lt.0.0.and.lat(ni,j).ge.0.0))
     &        lfndpoles=.true.
        enddo
c first pass: assume no dateline or pole boundaries
        rnorth = -90.
        south  = +90.
        west = +1000.
        east = -1000.
        do i = 1,ni
           do j=1,nj
               rnorth = max(rnorth,lat(i,j))
               south  = min(south ,lat(i,j))
               east   = max(east  ,lon(i,j))
               west   = min(west  ,lon(i,j))
           enddo
        enddo
c
c do we have a boundary to deal with?
        
        if(lfnddateline)then      !dateline: want west the min positive lon
           west = +1000.          !          want east the max negative lon
           east = -1000.
           do j=1,nj
              if(lon(1,j).gt.0.0)
     &        west = min(west,lon(1,j))
              if(lon(ni,j).lt.0.0)
     &        east = max(east,lon(ni,j))
           enddo
        elseif(lfndgrenwich)then  !Greenwich: we are set to go:
       
        endif

        if(lfndpolen)then
           rnorth=+90.
           do i=1,ni
              rnorth=min(rnorth,lat(i,nj))
           enddo
        endif

        if(lfndpoles)then
           south=-90.
           do i=1,ni
              south=max(south,lat(i,1))
           enddo
        endif
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

c       rnorth = min(rnorth, +90.)
c       south  = max(south , -90.)
c       east   = min(east  ,+180.)
c       west   = max(west  ,-180.)

        return
        end

