
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
            rnorth = max(rnorth,lat(i,1),lat(i,nj))
            south  = min(south ,lat(i,1),lat(i,nj))
            east   = max(east  ,lon(i,1),lon(i,nj))
            west   = min(west  ,lon(i,1),lon(i,nj))
        enddo ! i

        do j = 1,nj
            rnorth = max(rnorth,lat(1,j),lat(ni,j))
            south  = min(south ,lat(1,j),lat(ni,j))
            east   = max(east  ,lon(1,j),lon(ni,j))
            west   = min(west  ,lon(1,j),lon(ni,j))
        enddo ! j

        rnorth = rnorth + r_buffer
        south  = south  - r_buffer
        east   = east   + r_buffer
        west   = west   - r_buffer

        write(6,101)rnorth,south,east,west
101     format(1x,' Box around LAPS grid - NSEW ',4f9.2)

        return
        end

