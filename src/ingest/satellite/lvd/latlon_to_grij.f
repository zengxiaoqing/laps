
        subroutine latlon_to_grij(lat,lon,nx_l,ny_l,
     1                            lats,lons,gri,grj,nx_s,ny_s)

!       Determine i/j satellite coordinates for each LAPS grid point

        real lat_l(nx_l,ny_l)
        real lon_l(nx_l,ny_l)
        real lat_s(nx_s,ny_s)
        real lon_s(nx_s,ny_s)

        real ri_laps(nx_s,ny_s)
        real rj_laps(nx_s,ny_s)

        real a(0:1,0:1)
        real ainv(0:1,0:1)
        real b(0:1)
        real bsol(0:1)

!       Start with i/j LAPS coordinate for each satellite grid point

        do is = 1,nx_s
        do js = 1,ny_s
            call latlon_to_rlapsgrid(lats(is,js),lons(is,js),lat,lon
     1                              ,nx_l,ny_l
     1                              ,ri_laps(is,js),rj_laps(is,js)
     1                              ,istatus)
        enddo ! jl
        enddo ! il

        n = 2
        m = 2

        iwrite = 1

        do il = 1,nx_l
        do jl = 1,ny_l

            is_start = 2
            is_end = nx_s
            js_start = 2
            js_end = ny_s

            ri_s = 5.
            rj_s = 5.

            is = int(ri_s)
            js = int(rj_s)

            do iter = 1,1

!               Determine if LAPS grid point is within satellite box
                dlapsi_dx = ri_laps(is+1,js) - ri_laps(is,js)
                dlapsi_dy = ri_laps(is,js+1) - ri_laps(is,js)
                dlapsj_dx = rj_laps(is+1,js) - rj_laps(is,js)
                dlapsj_dy = rj_laps(is,js+1) - rj_laps(is,js)
                a(0,0) = dlapsi_dx
                a(0,1) = dlapsi_dy
                a(1,0) = dlapsj_dx
                a(1,1) = dlapsj_dy
                b(0) = il
                b(1) = jl
                ainv = a
                bsol = b
                call gaussj(ainv,n,n,bsol,m,m,ierr)
                if(iwrite .eq. 1)then
                    write(6,*)'a = ',a
                    write(6,*)'ainv = ',ainv
                    write(6,*)'b = ',b
                    write(6,*)'bsol = ',bsol
                endif
                iwrite = 0
            enddo ! iter

        enddo ! jl
        enddo ! il

        return
        end
