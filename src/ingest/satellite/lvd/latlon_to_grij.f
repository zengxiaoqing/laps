
        subroutine latlon_to_grij(lat_l,lon_l,nx_l,ny_l,
     1                            lat_s,lon_s,gri,grj,nx_s,ny_s,istatus)

        ANGDIF(XX,Y)=MOD(XX-Y+540.,360.)-180.

!       Determine i/j satellite coordinates for each LAPS grid point

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

        real a(0:1,0:1)
        real ainv(0:1,0:1)
        real b(0:1)
        real bsol(0:1)
        real x(0:1)

        character*1 c_loop

        real bi_coeff(2,2)

        itstatus=ishow_timer()

        write(6,*)' Subroutine latlon_to_grij...'

        istatus = 1

        call get_r_missing_data(r_missing_data,istatus)

        write(6,*)' lat_s range ',minval(lat_s),maxval(lat_s)
        write(6,*)' lat_l range ',minval(lat_l),maxval(lat_l)

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
        do is = 1,nx_s
        do js = 1,ny_s
            if(lat_s(is,js) .ne. r_missing_data .AND.
     1         lon_s(is,js) .ne. r_missing_data       )then
                call latlon_to_rlapsgrid(lat_s(is,js),lon_s(is,js)
     1                                  ,lat_l,lon_l
     1                                  ,nx_l,ny_l
     1                                  ,rilaps_s(is,js),rjlaps_s(is,js)
     1                                  ,istatus)
                slatmin = min(slatmin,lat_s(is,js))
                slatmax = max(slatmax,lat_s(is,js))

                slondiff = angdif(lon_s(is,js),sloncen)

                slonmin = min(slonmin,slondiff)
                slonmax = max(slonmax,slondiff)
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

!       Create buffer within sat lat/lon range
        sat_ll_buf = 0.2
        slatmin = slatmin + sat_ll_buf
        slatmax = slatmax - sat_ll_buf
        slonmin = slonmin + sat_ll_buf
        slonmax = slonmax - sat_ll_buf

        n = 2
        m = 2

        iwrite = 1
        iclo_call = 0
        iclo_loop = 0
        iclo_miss = 0
        iclo_conv = 0
        icall_gaussj = 0

        itstatus=ishow_timer()

!       Obtain approximate satellite geometry info
        write(6,*)' Calling satgeom...'
        range_m = 42155680.00
        call satgeom(i4time,lat_l,lon_l,nx_l,ny_l
     1  ,sublat_l,sublon_l,range_m,r_missing_data,Phase_angle_d
     1  ,Specular_ref_angle_d,emission_angle_d,azimuth_d,istatus)

        itstatus=ishow_timer()

        nxl2 = nx_l / 2  

!       Loop through LAPS grid points
        do il = 1,nx_l

          if(il .eq. (il/100) * 100)then
            write(6,*)'iclo_call = ',iclo_call
            write(6,*)'iclo_loop = ',iclo_loop
            write(6,*)'iclo_miss = ',iclo_miss
            write(6,*)'iclo_conv = ',iclo_conv
            write(6,*)'icall_gaussj = ',icall_gaussj
            itstatus=ishow_timer()
          endif

          do jl = 1,ny_l

            angdist = sqrt(angdif(lon_l(il,jl),sloncen)**2 +
     1                           (lat_l(il,jl)-slatcen)**2  )

            if(il .eq. il/100*100 .and. jl .eq. jl/100*100)then
                idebug = 1
            else
                idebug = 0
            endif

            if(emission_angle_d(il,jl) .gt. 360.)then
                write(6,*)
     1              ' ERROR in latlon_to_grij.f: large emission angle'
     1             ,emission_angle_d(il,jl)
                istatus = 0
                return
            endif

            if(idebug .eq. 1)then
                write(6,*)'Debugging - reached ',il,jl
                write(6,51)lat_l(il,jl),lon_l(il,jl)
     1                    ,angdist,emission_angle_d(il,jl)
51              format(' lat/lon/angdist/emission = ',2f9.3,2x,2f9.3)  
            endif

!           Reset when starting a new column
            if(jl .eq. 1)then
                bnorm = 999.
            endif

!           Test whether in satellite sector / geometry area
            if(lat_l(il,jl) .ge. slatmin .AND. 
     1         lat_l(il,jl) .le. slatmax .AND.
     1         il .gt. 1 .and. il .lt. nx_l .AND.
     1         angdif(lon_l(il,jl),sloncen) .ge. slonmin .AND.       
     1         angdif(lon_l(il,jl),sloncen) .le. slonmax .AND.
     1                 emission_angle_d(il,jl) .ge. 13.
!    1                 angdist                 .le. 70.
     1                                                           )then

              do iter = 1,5

                if(bnorm .gt. 1.1)then
                    
!                   Determine satellite grid point having the closest lat/lon
!                   to the LAPS grid point

!                   Apply initial check prior to routine to reduce calling
!                   Check if we've moved far enough since the last full calculation
                    delt_pt = sqrt((rlat_last-lat_l(il,jl))**2 
     1                           + (rlon_last-lon_l(il,jl))**2)
                    if((dist_min_last - delt_pt) .gt. 1.0)then
                       i_loop = 0
                       continue ! return
                    else 
                       if(.false.)then
                          call get_closest_point(
     1                             lat_l(il,jl),lon_l(il,jl)
     1                            ,r_missing_data           
     1                            ,lat_s,lon_s,nx_s,ny_s 
     1                            ,iclo_loop,i_loop                  ! I/O
     1                            ,rlat_last,rlon_last,dist_min_last ! I/O
     1                            ,iclo,jclo,dist_min)               ! O

                       else ! inlined version for speed
!                         Check if we've moved far enough since the last
!                         full calculation
                          rlat = lat_l(il,jl)
                          rlon = lon_l(il,jl)
                          delt_pt = sqrt((rlat_last-rlat)**2 
     1                                 + (rlon_last-rlon)**2)
                          if((dist_min_last - delt_pt) .gt. 1.0)then
                            i_loop = 0
                            goto 100
                          endif

                          iclo = 0
                          jclo = 0

                          dist_min = 999.

                          scale_lon = cosd(rlat)

                          do i = 1,nx_s
                          do j = 1,ny_s
                            if(lat_s(i,j) .ne. r_missing_data .AND.
     1                         lon_s(i,j) .ne. r_missing_data      )then
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
100                       continue
                       endif ! end of inlined code
                    endif

                    if(i_loop .eq. 1)then
                        c_loop = '*'
                    else
                        c_loop = ' '
                    endif

                    iclo_call = iclo_call + 1
                    if(iclo_call .le. 500 .or. 
     1                 dist_min .gt. dist_min_thr)       then
                        if(bnorm .ne. 20000.)then
                            write(6,101)iclo_call,iclo_loop,c_loop
     1                                     ,il,jl,iclo,jclo
     1                                     ,lat_l(il,jl),lon_l(il,jl)
     1                                     ,iter,dist_min,bnorm
101                         format(
     1    ' iclo_call/loop,il,jl,iclo,jclo,lat/lon,iter,dist_min,bnorm '       
     1                   ,2i6,1x,a1,1x,2i6,2x,2i6,2f10.4,i4,f11.4,f14.5)
                        else
                            write(6,102)iclo_call,iclo_loop,c_loop
     1                                     ,il,jl,iclo,jclo
     1                                     ,lat_l(il,jl),lon_l(il,jl)
     1                                     ,iter,dist_min,bnorm
102                         format(
     1    ' iclo_call/loop,il,jl,iclo,jclo,lat/lon,iter,dist_min,bnorm '       
     1                   ,2i6,1x,a1,1x,2i6,2x,2i6,2f10.4,i4,f11.4,f14.5
     1                          ,' no bilinear_laps result?')
                        endif
                    endif

                    if(dist_min .gt. dist_min_thr)then
                        gri(il,jl) = r_missing_data
                        grj(il,jl) = r_missing_data
                        iclo_miss = iclo_miss + 1
                        goto 900
                    endif

!                   Use closest point as first guess
                    ri_s = iclo
                    rj_s = jclo

                endif ! bnorm > 1.1

!               Interpolate to get laps i,j at current guessed sat i,j
                i1 = max(min(int(ri_s),nx_s-1),1); fi = ri_s - i1; i2=i1+1
                j1 = max(min(int(rj_s),ny_s-1),1); fj = rj_s - j1; j2=j1+1

                bi_coeff(1,1) = (1.-fi) * (1.-fj)
                bi_coeff(2,1) = fi      * (1.-fj)
                bi_coeff(1,2) = (1.-fi) *     fj 
                bi_coeff(2,2) = fi      *     fj 

                result1 = sum(bi_coeff(:,:) * rilaps_s(i1:i2,j1:j2))
                result2 = sum(bi_coeff(:,:) * rjlaps_s(i1:i2,j1:j2))

!               call bilinear_laps(ri_s,rj_s,nx_s,ny_s,rilaps_s,result1)
!               call bilinear_laps(ri_s,rj_s,nx_s,ny_s,rjlaps_s,result2)

!               call bilinear_interp_extrap(ri_s,rj_s,nx_s,ny_s
!    1                                     ,rilaps_s,result1,istat_bil)
!               call bilinear_interp_extrap(ri_s,rj_s,nx_s,ny_s
!    1                                     ,rjlaps_s,result2,istat_bil)

                if(result1 .eq. r_missing_data .OR. 
     1             result2 .eq. r_missing_data      )then
                    bnorm = 20000. ! flag value
                else
                    b(0) = il - result1

!                   Check for wrapping
                    if(b(0) .gt. (nxl2))then
                        b(0) = b(0) - nx_l
                    elseif(b(0) .lt. -(nxl2))then
                        b(0) = b(0) + nx_l
                    endif 

                    b(1) = jl - result2
                    bnorm = maxval(abs(b))
                endif

                is = int(ri_s)
                js = int(rj_s)

!               write(6,*)'is,js,ri_s,rj_s',is,js,ri_s,rj_s

                if(bnorm .le. epsilon)then ! convergence acheived
                    gri(il,jl) = ri_s             
                    grj(il,jl) = rj_s             
                    iclo_conv = iclo_conv + 1
                    if(iwrite .le. 60 .or. idebug .eq. 1)then
                        write(6,111)il,jl,iter,ri_s,rj_s
     1                             ,lat_l(il,jl),lon_l(il,jl)
111                     format(' Convergence acheived at ',3i6,2f10.4
     1                                                    ,4x,2f10.3)
                        if(idebug .eq. 1)write(6,*)
                        iwrite = iwrite + 1               
                    endif
                    goto 900
                elseif(bnorm .gt. 1000.)then
                    goto 800
                endif

!               Solve Ax = b for offsets of i/j satellite grid point
!               A matrix is partials of LAPS i,j with respect to satellite i,j
!               x matrix is the solution offset of satellite i,j
!               b matrix is desired LAPS i,j offset

!               Prevent partials from going beyond the array edge
                iis = max(min(is,nx_s-1),1)
                jjs = max(min(js,ny_s-1),1)

!               if(iis .le. 0)then
!                   write(6,*)'ERROR: iis,is,js,ri_s,rj_s'
!    1                               ,iis,is,js,ri_s,rj_s
!               endif

                dlapsi_dx = rilaps_s(iis+1,jjs) - rilaps_s(iis,jjs)
                dlapsi_dy = rilaps_s(iis,jjs+1) - rilaps_s(iis,jjs)

!               Test for wrapping
                if(dlapsi_dx .gt. nxl2)then
                    dlapsi_dx = dlapsi_dx - nx_l
                elseif(dlapsi_dx .lt. -nxl2)then
                    dlapsi_dx = dlapsi_dx + nx_l
                endif

                if(dlapsi_dy .gt. nxl2)then
                    dlapsi_dy = dlapsi_dy - nx_l
                elseif(dlapsi_dy .lt. -nxl2)then
                    dlapsi_dy = dlapsi_dy + nx_l
                endif

                dlapsj_dx = rjlaps_s(iis+1,jjs) - rjlaps_s(iis,jjs)
                dlapsj_dy = rjlaps_s(iis,jjs+1) - rjlaps_s(iis,jjs)

                a(0,0) = dlapsi_dx
                a(0,1) = dlapsi_dy
                a(1,0) = dlapsj_dx
                a(1,1) = dlapsj_dy

                ainv = a
                bsol = b

                if(.true.)then
                    call gaussj(ainv,n,n,bsol,m,m,ierr)
                else
                    aa = dlapsi_dx
                    bb = dlapsi_dy
                    cc = dlapsj_dx
                    dd = dlapsj_dy
                    arg1 = aa * dd - bb * cc
                    if(arg1 .eq. 0.)then
                        continue
                    else
                        arg2 = 1. / arg1
                    endif
                    aainv =  dd * arg2
                    bbinv = -bb * arg2
                    ccinv = -cc * arg2
                    ddinv =  aa * arg2

!                   Note that X = ainv times b
                    bsol(0) = aainv * b(0) + bbinv * b(1)
                    bsol(1) = ccinv * b(0) + ddinv * b(1)
                endif

                icall_gaussj = icall_gaussj + 1

                x = bsol

                ri_s = ri_s + x(0)
                rj_s = rj_s + x(1)

                if(iwrite .le. 10 .or. idebug .eq. 1)then
                    write(6,*)
                    write(6,*)'iter = ',iter
                    write(6,*)'a = ',a
                    write(6,*)'ainv = ',ainv
                    write(6,*)'b = ',b
                    write(6,*)'bsol/x = ',bsol
                    write(6,*)'bnorm = ',bnorm
                endif

              enddo ! iter

800           continue

!             Convergence was not achieved
              gri(il,jl) = r_missing_data
              grj(il,jl) = r_missing_data

!             if(iwrite .le. 20 .or. idebug .eq. 1)then
              if(.true.)then                                 
                if(bnorm .eq. 20000.)then
                    write(6,801)il,jl,bnorm,lat_l(il,jl),lon_l(il,jl)
801                 format('No bilinear_laps result at  ',2i5,f14.5
     1                    ,2f8.2)
                else
                    write(6,802)il,jl,bnorm,lat_l(il,jl),lon_l(il,jl) 
     1                         ,emission_angle_d(il,jl)
802                 format('Convergence not acheived at ',2i5,f14.5
     1                    ,2f8.2,' emission = ',f8.2)
                    if(abs(lon_l(il,jl)) .lt. 179.50)then
                        write(6,*)' ERROR not due to dateline'
                        istatus = 0
                        return
                    endif
                endif
                if(idebug .eq. 1)write(6,*)
                iwrite = iwrite + 1               
              endif

            else ! outside satellite sector
              dist_min = 999.
              if(idebug .eq. 1)then
                  write(6,*)'   Outside lat/lon/geom range'
     1                     ,lat_l(il,jl),lon_l(il,jl)
              endif
 
              gri(il,jl) = r_missing_data
              grj(il,jl) = r_missing_data

            endif ! in satellite sector          

900         continue

          enddo ! jl
        enddo ! il

        if(.true.)then ! extrapolate for wrapped domain edges
          do jl = 1,ny_l
            if(gri(2,jl) .ne. r_missing_data .and.
     1         gri(3,jl) .ne. r_missing_data )then
                gri(1,jl) = 2.*gri(2,jl) - gri(3,jl)
                grj(1,jl) = 2.*grj(2,jl) - grj(3,jl)
            endif
            if(gri(nx_l-1,jl) .ne. r_missing_data .and.
     1         gri(nx_l-2,jl) .ne. r_missing_data )then
                gri(nx_l,jl) = 2.*gri(nx_l-1,jl) - gri(nx_l-2,jl)
                grj(nx_l,jl) = 2.*grj(nx_l-1,jl) - grj(nx_l-2,jl)
            endif
          enddo ! jl
        endif

        write(6,*)'iclo_call = ',iclo_call
        write(6,*)'iclo_loop = ',iclo_loop
        write(6,*)'iclo_miss = ',iclo_miss
        write(6,*)'iclo_conv = ',iclo_conv
        write(6,*)'icall_gaussj = ',icall_gaussj

        return
        end
 
        subroutine get_closest_point(rlat,rlon,r_missing_data	
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
