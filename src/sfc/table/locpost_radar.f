
	subroutine locpost_radar(ni,nj,lat,lon,topo,ldf
     1                          ,lun,lun_out,istatus)

        character*200 static_dir,cfg_fname,line,radars_fn
        character*4 c4_name
        character*6 c6_chars
        character*7 c7_chars

        real lat(ni,nj),lon(ni,nj)
        real topo(ni,nj),ldf(ni,nj)

        write(6,*)' Subroutine locpost_radar...'

!       Get actual grid spacing valid at the gridpoint nearest the center
        icen = ni/2 + 1
        jcen = nj/2 + 1
        call get_grid_spacing_actual_xy(lat(icen,jcen),lon(icen,jcen)       
     1                        ,grid_spacing_actual_mx
     1                        ,grid_spacing_actual_my
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error return from get_grid_spacing_actual_xy'       
            stop
        endif

        grid_spacing_m = grid_spacing_actual_my


        call get_directory('static',static_dir,len_dir)

        cfg_fname = static_dir(1:len_dir)//'/NexradSite.cfg'
        call s_len(cfg_fname,len_fname)
        open(lun,file=cfg_fname(1:len_fname),status='old',err=998)

        radars_fn = static_dir(1:len_dir)//'/radarlist.dat'
        call s_len(radars_fn,len_rname)
        open(lun_out,file=radars_fn(1:len_rname),status='unknown'
     1                                          ,err=998)

        do ih = 1,13
            read(lun,*)
        enddo

        write(lun_out,*)'radar  perimeter  ocean      i      j'

        do il = 1,200
            read(lun,11,err=999,end=999)line
 11	    format(a)

            inum = 0
            do ic = 1,100
                if(line(ic:ic) .eq. '#')then
                    inum = inum + 1
                    if(inum .eq. 1)then
                        continue
                    elseif(inum .eq. 2)then
                        c4_name = line(ic+1:ic+4)
                    elseif(inum .eq. 3)then
                        continue
                    elseif(inum .eq. 4)then
                        c6_chars = line(ic+1:ic+6)
                        read(c6_chars,12)id,im,is
 12			format(3i2)
                        rlat = float(id) + float(im) / 60. 
     1                       + float(is)/3600.
                    elseif(inum .eq. 5)then
                        c7_chars = line(ic+1:ic+7)
                        read(c7_chars,13)id,im,is
 13			format(i3,i2,i2)
                        rlon = -float(id) + float(im) / 60. 
     1                        + float(is)/3600.
                    endif
                endif
            enddo ! ic

            call latlon_to_rlapsgrid(rlat,rlon,lat,lon
     1                              ,ni,nj,ri,rj
     1                              ,istatus)

            if(ri .le. 1.)then
                dist_outside_i = 1. - ri
            elseif(ri .gt. float(ni))then
                dist_outside_i = ri - float(ni)
            else
                dist_outside_i = 0.
            endif              

            if(rj .le. 1.)then
                dist_outside_j = 1. - rj
            elseif(rj .gt. float(nj))then
                dist_outside_j = rj - float(nj)
            else
                dist_outside_j = 0.
            endif              

            dist_outside = sqrt(dist_outside_i**2 + dist_outside_j**2)       

            write(6,*)line
            write(6,*)c4_name,rlat,rlon,ri,rj,dist_outside

            dist_outside_km = dist_outside * grid_spacing_m / 1000.

!           Find nearest ocean to radar site
            grid_dist_min = 9999.
            do ii = 1,ni
            do jj = 1,nj
                if(ldf(ii,jj) .lt. .01 .and. topo(ii,jj) .lt. 10.)then
                    grid_dist = sqrt((ri-ii)**2 + (rj-jj)**2) 
     1                          * (grid_spacing_m / 1000.)
                    grid_dist_min = min(grid_dist,grid_dist_min)
                endif
            enddo ! jj
            enddo ! ii
            dist_ocean_min = grid_dist_min 

            if(dist_outside_km .le. 100.)then
                write(lun_out,21)c4_name,dist_outside_km,dist_ocean_min       
     1                          ,nint(ri),nint(rj)
 21		format(1x,a4,f9.2,2x,f9.1,2i7)
            endif

        enddo ! il

        go to 999

 998	write(6,*)' Error in locpost_radar'
        istatus = 0

 999	continue

        istatus = 1

        close(lun)
        close(lun_out)

        return
        end
