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

        real*4 m ! Grid points per meter

        DATA scale/1./

        real*4 lat(ni,nj),lon(ni,nj)
        real*4 uanl(ni,nj),vanl(ni,nj)
        real*4 uanl_grid(ni,nj),vanl_grid(ni,nj)
        real*4 div(ni,nj)
        real*4 one(ni,nj)

        real*4 dum1(ni,nj)
        real*4 dum2(ni,nj)
        real*4 dum3(ni,nj)

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


        subroutine vorticity_abs(uanl,vanl,vort,lat,lon,ni,nj
     1                          ,l_grid_north,r_missing_data)

        real*4 lat(ni,nj),lon(ni,nj)
        real*4 uanl(ni,nj),vanl(ni,nj)
        real*4 coriolis(ni,nj)
        real*4 vort(ni,nj)
        real*4 div(ni,nj)
        real*4 dx(ni,nj), dy(ni,nj)

        call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)

        call vortdiv(uanl,vanl,vort,div,ni,nj,dx,dy)

        call get_coriolis_rotation(ni,nj,lat,coriolis)

        call add(vort,coriolis,vort,ni,nj)

        return
        end


        subroutine get_coriolis_rotation(nx,ny,lat,coriolis_rotation)
c
        include 'trigd.inc'
        implicit none

        integer nx,ny
        integer i,j

        real*4  lat(nx,ny)
        real*4  coriolis_rotation(nx,ny)

        real*4  omega_ear
        data    omega_ear/7.292e-5/

        do j=1,ny
        do i=1,nx
           coriolis_rotation(i,j)=2*omega_ear*sind(lat(i,j))
        enddo
        enddo

        return
        end

