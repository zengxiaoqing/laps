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


        subroutine cpt_pcp_cnc(ref_3d,temp_3d,cldpcp_type_3d  ! Input
     1                                  ,ni,nj,nk     ! Input
     1                                  ,pcp_cnc_3d   ! Output (kg/m^3)
     1                                  ,rai_cnc_3d   ! Output (rain)
     1                                  ,sno_cnc_3d   ! Output (snow)
     1                                  ,pic_cnc_3d)  ! Output (precip ice)

        real*4 temp_3d(ni,nj,nk)
        real*4 ref_3d(ni,nj,nk)
        integer cldpcp_type_3d(ni,nj,nk)

        real*4 pcp_cnc_3d(ni,nj,nk)
        real*4 rai_cnc_3d(ni,nj,nk)
        real*4 sno_cnc_3d(ni,nj,nk)
        real*4 pic_cnc_3d(ni,nj,nk)

        real*4 rate_3d(ni,nj,nk)

!       Get 3D precip rates one layer at a time with ZR routine
        do k = 1,nk
            call zr(ref_3d(1,1,k),ni,nj,rate_3d(1,1,k))
        enddo ! k

        do j = 1,nj
        do i = 1,ni
        do k = 1,nk
            pressure = pressure_of_level(k)
            ipcp_type = cldpcp_type_3d(i,j,k) / 16   ! Pull out precip type
            if(ipcp_type .ne. 0)then
                call cpt_fall_velocity(ipcp_type,pressure,temp_3d(i,j,k)
     1                                                  ,fall_velocity)
                call cpt_concentration(rate_3d(i,j,k),fall_velocity
     1                                          ,pcp_cnc_3d(i,j,k))

                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    rai_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    sno_cnc_3d(i,j,k) = 0.
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 2)then                       ! snow
                    rai_cnc_3d(i,j,k) = 0.
                    sno_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! IP or Hail
                    rai_cnc_3d(i,j,k) = 0.
                    sno_cnc_3d(i,j,k) = 0.
                    pic_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)

                endif ! ipcp_type

            else  ! ipcp_type = 0
                pcp_cnc_3d(i,j,k) = 0.
                rai_cnc_3d(i,j,k) = 0.
                sno_cnc_3d(i,j,k) = 0.
                pic_cnc_3d(i,j,k) = 0.

            endif ! ipcp_type .ne. 0

        enddo ! k
        enddo ! i
        enddo ! j

        end

        subroutine cpt_fall_velocity(ipcp_type,p,t,fall_velocity
     1                              ,ni,nj,nk)

        if(ipcp_type .eq. 1)then ! Rain
            fall_velocity = 5.0
        elseif(ipcp_type .eq. 2)then ! Snow
            fall_velocity = 1.0
        elseif(ipcp_type .eq. 3)then ! Freezing Rain
            fall_velocity = 5.0
        elseif(ipcp_type .eq. 4)then ! Sleet
            fall_velocity = 5.0
        elseif(ipcp_type .eq. 5)then ! Hail
            fall_velocity = 10.0
        endif

        if(p .ne. 101300. .or. t .ne. 273.15)then
            density_norm = (p / 101300.) * (273.15 / t)
            sqrt_density_norm = sqrt(density_norm)
        else
            sqrt_density_norm = 1.
        endif

        rate = rate * sqrt_density_norm

        return
        end


        subroutine cpt_concentration(rate,fall_velocity,concentration)

!       This gets m**3 of water per m**3 of volume (non-dimensional)
        concentration = rate / fall_velocity

        rho = 1e3   ! mass in kilograms of a cubic meter of water

        concentration = concentration * rho

        return

        end
