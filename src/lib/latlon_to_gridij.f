
        subroutine latlon_to_grij_new(lat_l,lon_l,nx_l,ny_l,
     1                                lat_s,lon_s,gri,grj,nx_s,ny_s,
     1                                istatus)

        ANGDIF(XX,Y)=MOD(XX-Y+540.,360.)-180.

!       Determine i/j satellite/mdl coordinates for each LAPS grid point

        real lat_l(nx_l,ny_l)
        real lon_l(nx_l,ny_l)
        real lat_s(nx_s,ny_s)
        real lon_s(nx_s,ny_s)

        real sublat_l(nx_l,ny_l)
        real sublon_l(nx_l,ny_l)

        real emission_angle_d(nx_l,ny_l)
        real phase_angle_d(nx_l,ny_l)
        real specular_ref_angle_d(nx_l,ny_l)
        real azimuth_d(nx_l,ny_l)

        real gri(nx_l,ny_l)
        real grj(nx_l,ny_l)

        real rilaps_s(nx_s,ny_s)
        real rjlaps_s(nx_s,ny_s)

        character*1 c_loop

        real bi_coeff(2,2)

        xlin(a,b,p,c,d,q) = (p-a*q/c) / (b-a*d/c)
        ylin(a,b,p,c,d,q) = (q-c*p/a) / (d-c*b/a)

        itstatus=ishow_timer()

        write(6,*)' Subroutine latlon_to_grij...'

        istatus = 1

        call get_r_missing_data(r_missing_data,istatus)

!       Initialize
        gri = r_missing_data
        grj = r_missing_data

        write(6,*)' lat_s range ',minval(lat_s),maxval(lat_s)
        write(6,*)' lat_l range ',minval(lat_l),maxval(lat_l)

        write(6,*)' lon_s range ',minval(lon_s),maxval(lon_s)
        write(6,*)' lon_l range ',minval(lon_l),maxval(lon_l)

        slatmax = -999.
        slatmin = +999.
        slonmax = -999.
        slonmin = +999.

        slatcen = lat_s(nx_s/2,ny_s/2)
        sloncen = lon_s(nx_s/2,ny_s/2)

        sublat_l = slatcen
        sublon_l = sloncen

        dist_min_thr = 0.5

        epsilon = .00040 ! degrees

!       Initialize other variables
        iclo = -1
        jclo = -1
        dist_min_last = -999.

!       Start with i/j LAPS coordinate for each satellite grid point
        write(6,*)' Fill ri/rj LAPS values for each sat grid point'

        do is = 1,nx_s
        do js = 1,ny_s

            if(is .eq. nx_s/2 .and. js .eq. ny_s/2)then
              iwrite = 1
            else
              iwrite = 0
            endif

            if(lat_s(is,js) .ne. r_missing_data .AND.
     1         lon_s(is,js) .ne. r_missing_data       )then
                call latlon_to_rlapsgrid(lat_s(is,js),lon_s(is,js)
     1                                  ,lat_l,lon_l
     1                                  ,nx_l,ny_l
     1                                  ,rilaps_s(is,js),rjlaps_s(is,js)
     1                                  ,istatus)
                if(iwrite .eq. 1)then
                  write(6,*)is,js,lat_s(is,js),lon_s(is,js),
     1                      rilaps_s(is,js),rjlaps_s(is,js)
                endif

                if(istatus .eq. 1)then
                   slatmin = min(slatmin,lat_s(is,js))
                   slatmax = max(slatmax,lat_s(is,js))

                   slondiff = angdif(lon_s(is,js),sloncen)

                   slonmin = min(slonmin,slondiff)
                   slonmax = max(slonmax,slondiff)
                else
                   rilaps_s(is,js) = r_missing_data
                   rjlaps_s(is,js) = r_missing_data
                endif
            else
                rilaps_s(is,js) = r_missing_data
                rjlaps_s(is,js) = r_missing_data
            endif
        enddo ! js
        enddo ! is

        write(6,*)' Sat lon center                  : ',sloncen
        write(6,*)' Sat lat range                   : ',slatmin,slatmax
        write(6,*)' Sat lon range (diff from center): ',slonmin,slonmax
        write(6,*)' Sat lon range                   : ',slonmin+sloncen
     1                                                 ,slonmax+sloncen

        itstatus=ishow_timer()

!       Create buffer within sat lat/lon range
        sat_ll_buf = 0.2
        slatmin = slatmin + sat_ll_buf
        slatmax = slatmax - sat_ll_buf
        slonmin = slonmin + sat_ll_buf
        slonmax = slonmax - sat_ll_buf

!       New solution method
        do is = 1,nx_s-1
        do js = 1,ny_s-1

          if(is .eq. nx_s/2 .and. js .eq. ny_s/2)then
            iwrite = 1
          else
            iwrite = 0
          endif

!         Obtain derivatives
          dlidsi = rilaps_s(is+1,js) - rilaps_s(is,js)
          dlidsj = rilaps_s(is,js+1) - rilaps_s(is,js)
          dljdsi = rjlaps_s(is+1,js) - rjlaps_s(is,js)
          dljdsj = rjlaps_s(is,js+1) - rjlaps_s(is,js)

!         Locate model grid points inside the satellite parallelogram
          rilmn = min(rilaps_s(is,js),rilaps_s(is+1,js)
     1               ,rilaps_s(is,js+1),rilaps_s(is+1,js+1))

          rilmx = max(rilaps_s(is,js),rilaps_s(is+1,js)
     1               ,rilaps_s(is,js+1),rilaps_s(is+1,js+1))

          rjlmn = min(rjlaps_s(is,js),rjlaps_s(is+1,js)
     1               ,rjlaps_s(is,js+1),rjlaps_s(is+1,js+1))

          rjlmx = max(rjlaps_s(is,js),rjlaps_s(is+1,js)
     1               ,rjlaps_s(is,js+1),rjlaps_s(is+1,js+1))

          if(rilmn .ne. r_missing_data .and.
     1       rilmx .ne. r_missing_data .and.
     1       rjlmn .ne. r_missing_data .and.
     1       rjlmx .ne. r_missing_data       )then

            if(iwrite .eq. 1)then
               write(6,*)is,js,rilaps_s(is,js),rjlaps_s(is,js),
     1                             int(rilmn),int(rilmx),
     1                             int(rjlmn),int(rjlmx)  
            endif

!           Loop through model grid points inside the satellite parallelogram
            do il = int(rilmn)-1,int(rilmx)+1
            do jl = int(rjlmn)-1,int(rjlmx)+1

              if(il .ge. 1 .and. il .le. nx_l .and.
     1           jl .ge. 1 .and. jl .le. ny_l       )then

!               Set up system of linear equations
!               If x,y = 0 we are at the lower left of the satellite parallelogram
!               ax + by = p
!               cx + dy = q

                a = dlidsi 
                b = dlidsj
                c = dljdsi
                d = dljdsj

                p = float(il)         ! model grid point of interest
                q = float(jl)         ! model grid point of interest

                if(il .eq. 170 .and. jl .eq. 272)iwrite = 1

                if(a .ne. 0. .and. c .ne. 0.)then ! well conditioned
                    x = xlin(a,b,p,c,d,q) ! sat fractional grid point (fis)
                    y = ylin(a,b,p,c,d,q) ! sat fractional grid point (fjs)

                    gri(il,jl) = float(is) + x
                    grj(il,jl) = float(js) + y

                    dist = sqrt(x**2 + y**2)

                    if(iwrite .eq. 1)then
                      write(6,11)is,js,il,jl,gri(il,jl),grj(il,jl),dist
11                    format(' solved ',4i6,3f9.4)                   
                    endif

                endif

              endif ! inside laps grid

            enddo ! jl
            enddo ! il

          endif

        enddo ! js
        enddo ! is

        icount_good = 0
        icount_miss = 0

        do il = 1,nx_l
        do jl = 1,ny_l
          if(gri(il,jl) .ne. r_missing_data .and.
     1       grj(il,jl) .ne. r_missing_data      )then
            icount_good = icount_good + 1
          else
            icount_miss = icount_miss + 1
            write(6,*)' missing at ',il,jl
          endif
        enddo ! jl
        enddo ! il

        write(6,*)' Number of good/miss points is ',icount_good
     1                                             ,icount_miss

        itstatus=ishow_timer()
        
        return
        end
 
        subroutine get_closest_point_multiscale
     1                              (rlat,rlon,r_missing_data	
     1                              ,lat_s,lon_s,nx_s,ny_s 
     1                              ,iclo_loop,i_loop
     1                              ,rlat_last,rlon_last,dist_min_last ! I/O
     1                              ,iclo,jclo,dist_min)                

        include 'trigd.inc'

        real lat_s(nx_s,ny_s)
        real lon_s(nx_s,ny_s)

!       Determine satellite grid point having the closest lat/lon
!       to the LAPS grid point

!       Check if we've moved far enough since the last full calculation
        delt_pt = sqrt((rlat_last-rlat)**2 + (rlon_last-rlon)**2)
        if((dist_min_last - delt_pt) .gt. 1.0)then
            i_loop = 0
            return
        endif

        iclo = 0
        jclo = 0

        dist_min = 999.

        scale_lon = cosd(rlat)

        do i = 1,nx_s
        do j = 1,ny_s
          if(lat_s(i,j) .ne. r_missing_data .AND.
     1       lon_s(i,j) .ne. r_missing_data       )then
            delta_lat = abs(lat_s(i,j) - rlat)
            delta_lon = abs(lon_s(i,j) - rlon) * scale_lon
            dist_deg = sqrt(delta_lat**2 + delta_lon**2)
            if(dist_deg .lt. dist_min)then
                dist_min = dist_deg
                iclo = i
                jclo = j
            endif
          endif
        enddo ! j
        enddo ! i

        iclo_loop = iclo_loop + 1
        i_loop = 1

        rlat_last = rlat
        rlon_last = rlon
        dist_min_last = dist_min

        return
        end 

