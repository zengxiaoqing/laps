cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis

        subroutine divergence(uanl,vanl,div,lat,lon,ni,nj
     1  ,one,dum1,dum2,dum3,uanl_grid,vanl_grid,r_missing_data)

!      ~90            Steve Albers  Original Version
!       97-Aug-17     Ken Dritz     Added r_missing_data as dummy argument
!       97-Aug-17     Ken Dritz     Removed include of lapsparms.for
!       97-Oct        Steve Albers  Add lon in call to fflxc.

        include 'trigd.inc'

        real m ! Grid points per meter

        real one(ni,nj)

        DATA scale/1./

        real*4 lat(ni,nj),lon(ni,nj)
        real*4 uanl(ni,nj),vanl(ni,nj)
        real*4 uanl_grid(ni,nj),vanl_grid(ni,nj)
        real*4 div(ni,nj)

        PHI0 = standard_latitude

        grid_spacing_m = sqrt(
     1                 (  lat(1,2) - lat(1,1)                  )**2
     1               + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2
     1                                  )    * 111317. ! Grid spacing m
        m = 1.0 / grid_spacing_m

        do j = 1,nj
        do i = 1,ni
            uanl(i,j) = 0.
            vanl(i,j) = 500.0
            one(i,j) = 1.0
            call uvtrue_to_uvgrid(uanl(i,j),vanl(i,j)
     1          ,uanl_grid(i,j),vanl_grid(i,j),lon(i,j))
        enddo ! i
        enddo ! j

        call FFLXC(ni,nj,M,PHI0,SCALE
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
