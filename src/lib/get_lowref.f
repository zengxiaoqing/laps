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
        subroutine get_low_ref(z_3d_in,pres_sfc_pa,ni,nj,nk
     1                        ,radar_2d_out)

!       Steve Albers            1990

        real z_3d_in(ni,nj,nk)
        real radar_2d_out(ni,nj)
        real pres_sfc_pa(ni,nj)

        write(6,*)' Converting from 3D Z field to Low Level Z field'

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_low_ref, STOP'
            stop
        endif

        n_low_pts = 0

        do j = 1,nj
        do i = 1,ni

!           k_topo = max(int(height_to_zcoord(topo(i,j),istatus)),1)

            k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
            radar_2d_out(i,j) = z_3d_in(i,j,k_topo+1)
            if(radar_2d_out(i,j) .ne. ref_base)then
                n_low_pts = n_low_pts + 1
            endif

        enddo
        enddo

        write(6,*)' n_low_pts  = ',n_low_pts

        return
        end
