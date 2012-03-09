
	subroutine locpost_radar(ni,nj,lat,lon,topo,ldf     
     1                          ,lun,lun_out,istatus)

        character*200 static_dir,cfg_fname,line,radars_fn,out_fn
        character*100 c_remap
        character*800 c_wideband
        character*4 c4_name, c4_name_a(200), c4_name_lc
        character*7 c_vxx
        character*6 c6_chars
        character*7 c7_chars

        real lat(ni,nj),lon(ni,nj)
        real topo(ni,nj),ldf(ni,nj)

        write(6,*)' Subroutine locpost_radar...'

        ie = 0

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

!       Open NexradSite.cfg for reading
!       cfg_fname = static_dir(1:len_dir)//'/NexradSite.cfg'
        cfg_fname = static_dir(1:len_dir)//'/NexradSite.sorted'
        call s_len(cfg_fname,len_fname)
        open(lun,file=cfg_fname(1:len_fname),status='old',err=998)

!       Open radarlist.dat for output
        radars_fn = static_dir(1:len_dir)//'/radarlist.dat'
        call s_len(radars_fn,len_rname)
        open(lun_out,file=radars_fn(1:len_rname),status='unknown'
     1                                          ,err=998)

!       Open widebandlist_radars.txt for output
        lun_wideband = 41
        out_fn = static_dir(1:len_dir)//'/widebandlist_radars.txt'
        call s_len(out_fn,len_name)
        open(lun_wideband,file=out_fn(1:len_name),status='unknown'
     1                                           ,err=998)

        do ih = 1,13
            read(lun,*)
        enddo

        write(lun_out,*)'radar  perimeter  ocean      i      j'

        icount = 0

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
                        rlon = -float(id) - float(im) / 60. ! Terms are all negative for west longitude
     1                        - float(is)/3600.
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

            dist_outside_km = dist_outside * grid_spacing_m / 1000.

            write(6,*)line
            write(6,*)c4_name,rlat,rlon,ri,rj,dist_outside_km

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

            if(dist_outside_km .le. 200.)then ! apply distance threshold
                icount = icount + 1
                write(lun_out,21)c4_name,dist_outside_km,dist_ocean_min       
     1                          ,nint(ri),nint(rj)
 21		format(1x,a4,f9.2,2x,f9.1,2i7)

!               Write to string for widebandlist
                is = (icount-1) * 5 + 1
                ie = is + 4
                c_wideband(is:ie) = c4_name//' '
                write(6,*)' c_wideband section: ',c_wideband(is:ie)

!               Save string into array for later use
                c4_name_a(icount) = c4_name

            endif

        enddo ! il

        go to 999

 998	write(6,*)' Error in locpost_radar'
        istatus = 0

 999	continue

        if(ie .eq. 0)goto 9999

!       Open remap_radars.nl for output
        lun_remap = 51
        out_fn = static_dir(1:len_dir)//'/remap_radars.nl'           
        call s_len(out_fn,len_name)
        open(lun_remap,file=out_fn(1:len_name),status='unknown'
     1                                        ,err=998)

!       Write widebandlist to file
        write(6,*)                               
        write(6,*)' Full wideband string...'
        write(6,101)c_wideband(1:ie-1)
        write(lun_wideband,101)c_wideband(1:ie-1)
101     format(a)

!       Write to remap_radars.nl (n_radars_remap)    
        if(icount .le. 9)then
            write(lun_remap,102)icount
102         format(' n_radars_remap=',i1,',')
        elseif(icount .le. 99)then
            write(lun_remap,103)icount
103         format(' n_radars_remap=',i2,',')
        elseif(icount .le. 999)then
            write(lun_remap,104)icount
104         format(' n_radars_remap=',i3,',')
        endif

!       Write to remap_radars.nl (upper case section)
        write(lun_remap,*)        
        write(lun_remap,*)'path_to_radar_uc_a='
        do i = 1,icount ! upper case section
            write(c_remap,111)c4_name_a(i)
111         format("'/widebandroot/",a,"/netcdf',")
            write(lun_remap,*)trim(c_remap)
        enddo ! i

!       Write to remap_radars.nl (lower case section)
        write(lun_remap,*)
        write(lun_remap,*)'path_to_radar_lc_a='
        do i = 1,icount ! upper case section
            call downcase(c4_name_a(i),c4_name_lc)
            write(c_remap,111)c4_name_lc  
            write(lun_remap,*)trim(c_remap)
        enddo ! i

!       Write to remap_radars.nl (vxx section)
        write(lun_remap,*)
        write(lun_remap,*)'laps_radar_ext_a='
        do i = 1,icount ! upper case section
            call downcase(c4_name_a(i),c4_name_lc)
            if(icount .lt. 100)then
                write(c_vxx,121)i              
121             format("'v",i2.2,"',")
            else
                write(c_vxx,122)               
122             format("'v",i3.3,"',")
            endif
            write(lun_remap,*)trim(c_vxx)         
        enddo ! i

        write(lun_remap,*)'/'

        istatus = 1

 9999	close(lun)
        close(lun_out)
        close(lun_wideband)
        close(lun_remap)

        return
        end
