cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

        subroutine divergence(uanl,vanl,div,lat,lon,ni,nj
     1                       ,l_grid_north,r_missing_data)

!      ~90            Steve Albers  Original Version
!       97-Aug-17     Ken Dritz     Added r_missing_data as dummy argument
!       97-Aug-17     Ken Dritz     Removed include of lapsparms.for
!       97-Oct        Steve Albers  Add lon in call to fflxc.

        include 'trigd.inc'

        real m ! Grid points per meter

        DATA scale/1./

        real lat(ni,nj),lon(ni,nj)
        real uanl(ni,nj),vanl(ni,nj)
        real uanl_grid(ni,nj),vanl_grid(ni,nj)
        real div(ni,nj)
        real one(ni,nj)

        real dum1(ni,nj)
        real dum2(ni,nj)
        real dum3(ni,nj)

        character*6 c6_maproj

        logical l_grid_north

!       grid_spacing_m = sqrt(
!    1                 (  lat(1,2) - lat(1,1)                  )**2
!    1               + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2
!    1                                  )    * 111317. ! Grid spacing m

        call get_standard_latitudes(slat1,slat2,istatus)
        if(istatus .ne. 1)then
            return
        endif

        call get_grid_spacing(grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            return
        endif

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            return
        endif

        if(c6_maproj .eq. 'plrstr')then
            call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
        else
            grid_spacing_proj_m = grid_spacing_m
        endif

        write(6,*)' Grid spacings (m) = ',grid_spacing_m
     1                                   ,grid_spacing_proj_m

        m = 1.0 / grid_spacing_proj_m

        do j = 1,nj
        do i = 1,ni
!           uanl(i,j) = 0.
!           vanl(i,j) = 500.0
            one(i,j) = 1.0

            if(l_grid_north)then
                uanl_grid(i,j) = uanl(i,j)
                vanl_grid(i,j) = vanl(i,j)
            else
                call uvtrue_to_uvgrid(uanl(i,j),vanl(i,j)
     1                               ,uanl_grid(i,j),vanl_grid(i,j)
     1                               ,lon(i,j))
            endif

        enddo ! i
        enddo ! j

        call FFLXC(ni,nj,M,SCALE
     1  ,uanl_grid,vanl_grid,one,div,lat,lon
     1  ,dum1,dum2,dum3,r_missing_data)

        do j = 1,nj
        do i = 1,ni

            if(abs(div(i,j)) .gt. 1e10)then
                div(i,j) = 0.
            else
                div(i,j) = -div(i,j)
            endif

        enddo ! j
        enddo ! i

        return
        end


        subroutine vorticity_abs(uanl,vanl,vort,lat,lon,ni,nj,dx,dy
     1                          ,l_grid_north,r_missing_data)

        real lat(ni,nj),lon(ni,nj)                           ! I
        real uanl(ni,nj),vanl(ni,nj)                         ! I
        real coriolis(ni,nj)                                 ! L
        real vort(ni,nj)                                     ! O
        real div(ni,nj)                                      ! L
        real dx(ni,nj), dy(ni,nj)                            ! O

        logical l_grid_north                                 ! I

        data init/0/
        save init

        if(init .eq. 0)then ! call this only once to save time
            write(6,*)'     call get_grid_spacing_array'
            call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)
            init = 1
        endif

        write(6,*)'     call vortdiv'
        call vortdiv(uanl,vanl,vort,div,ni,nj,dx,dy)

        write(6,*)'     call get_coriolis_rotation'
        call get_coriolis_rotation(ni,nj,lat,coriolis)

!       call add(vort,coriolis,vort,ni,nj)
        vort = vort + coriolis 

        return
        end


        subroutine get_coriolis_rotation(nx,ny,lat,coriolis_rotation)
c
        include 'trigd.inc'
        implicit none

        integer nx,ny
        integer i,j

        real  lat(nx,ny)
        real  coriolis_rotation(nx,ny)

        real  omega_ear
        data    omega_ear/7.292e-5/

        do j=1,ny
        do i=1,nx
           coriolis_rotation(i,j)=2*omega_ear*sind(lat(i,j))
        enddo
        enddo

        return
        end



        subroutine calc_potvort(i4time,uanl,vanl,temp_3d,potvort,lat,lon       
     1                         ,ni,nj,nk,nkuv,k_in,l_grid_north,dx,dy
     1                         ,r_missing_data,istatus)       

        include 'constants.inc'

        real lat(ni,nj),lon(ni,nj)                           ! I
        real uanl(ni,nj,nkuv),vanl(ni,nj,nkuv)               ! I
        real temp_3d(ni,nj,nk)                               ! I
        real vort_2d(ni,nj)                                  ! L
        real dx(ni,nj)                                       ! O
        real dy(ni,nj)                                       ! O
        real theta(ni,nj,nk)                                 ! L
        real pres_3d(ni,nj,nk)                               ! L
        real dtheta_dp(ni,nj,nk)                             ! L
        real potvort(ni,nj,nk)                               ! O

        logical l_grid_north                                 ! I

        write(6,*)' Subroutine calc_potvort'

!       Obtain 3D pressure field
        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)return

!       Calculate theta
        do k = 1,nk
        do i = 1,ni
        do j = 1,nj
            theta(i,j,k) = o(temp_3d(i,j,k),pres_3d(i,j,k))
        enddo ! j
        enddo ! i
        enddo ! k

!       Calculate dtheta / dp
        do k = 2,nk-1
            kp1 = k+1
            km1 = k-1
            do i = 1,ni
            do j = 1,nj
                dtheta = theta(i,j,kp1) - theta(i,j,km1)              
                dp     = pres_3d(i,j,kp1) - pres_3d(i,j,km1)
                dtheta_dp(i,j,k) = dtheta / dp
            enddo ! j
            enddo ! i
        enddo ! k
        dtheta_dp(:,:,1)  = dtheta_dp(:,:,2)
        dtheta_dp(:,:,nk) = dtheta_dp(:,:,nk-1)


!       Calculate absolute vorticity
        if(k_in .gt. 0)then ! 2D
            kstart = k_in
            kend = k_in
        else                ! 3D where k_in = 0
            kstart = 1
            kend = nk
        endif

        do k = kstart,kend
            if(k_in .gt. 0)then ! 2D
                kuv = 1
            else                ! 3D where k_in = 0
                kuv = k
            endif

            write(6,*)' calling vorticity_abs for kuv = ',kuv
            call vorticity_abs(uanl(1,1,kuv),vanl(1,1,kuv),vort_2d
     1                   ,lat,lon,ni,nj,dx,dy
     1                   ,l_grid_north,r_missing_data)

!           See http://www-das.uwyo.edu/~geerts/cwx/notes/chap12/pot_vort.html

            if(k_in .gt. 0)then     ! just returning 2d field (for efficiency)
                write(6,*)'     calculating 2D potvort - with debugging'       
                do i = 1,ni
                do j = 1,nj
                    potvort(i,j,1) = vort_2d(i,j) * dtheta_dp(i,j,k) 
     1                                            * grav       
                    if(i .eq. ni/2)then ! debug
                        write(6,*)i,j,vort_2d(i,j)
     1                           ,uanl(i,j,1),vanl(i,j,1)
     1                           ,dtheta_dp(i,j,k)
     1                           ,potvort(i,j,1)
                    endif
                enddo ! j
                enddo ! i

            elseif(k_in .eq. 0)then ! return 3D field
                write(6,*)'     calculating 3D potvort'
                potvort(:,:,k) = vort_2d(:,:) * dtheta_dp(:,:,k) 
     1                                        * grav       
            endif
        enddo ! k

        return
        end
