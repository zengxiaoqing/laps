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

        subroutine helicity_laps(uanl,vanl,ustorm,vstorm
     1                          ,heights_3d,topo
     1                          ,imax,jmax,kmax,helicity,istatus)       

        real*4 ustorm(imax,jmax),vstorm(imax,jmax)
        real*4 usfc(imax,jmax),vsfc(imax,jmax)
        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)
        real*4 helicity(imax,jmax)
        real*4 topo(imax,jmax)

!       1998 Steve Albers - Overhauled

        icount_write = 0

        write(6,*)' Computing Helicity for 0-3 km AGL'

        do j = 1,jmax
        do i = 1,imax

            area_sum = 0.

!           Layer is 0-3 km AGL, denoted from "sfc" to "top"
            rksfc = height_to_zcoord2(topo(i,j)      ,heights_3d
     1                               ,imax,jmax,kmax,i,j,istatus)
            if(istatus .ne. 1)return

            rktop = height_to_zcoord2(topo(i,j)+3000.,heights_3d
     1                               ,imax,jmax,kmax,i,j,istatus)
            if(istatus .ne. 1)return

!           Get storm relative wind at the sfc, using interpolated sfc wind
            ksfc  = int(rksfc)
            frack = rksfc - ksfc
            u_sfc = uanl(i,j,ksfc)   * (1.-frack) 
     1            + uanl(i,j,ksfc+1) * frack       
            v_sfc = vanl(i,j,ksfc)   * (1.-frack) 
     1            + vanl(i,j,ksfc+1) * frack

            u_rel_l = u_sfc - ustorm(i,j)
            v_rel_l = v_sfc - vstorm(i,j)

            klow  = int(rksfc) + 1          ! 1st level above the sfc
            khigh = int(rktop) + 1          ! 1st level above top of layer

            do k = klow,khigh

                if(k .lt. khigh)then
!                   Get storm relative wind at this level
                    u_rel_u = uanl(i,j,k) - ustorm(i,j)
                    v_rel_u = vanl(i,j,k) - vstorm(i,j)

                else ! k = khigh, use 3km agl values instead of this laps level
                    frack = rktop - int(rktop)
                    u_anl = uanl(i,j,k-1) * (1.-frack) 
     1                    + uanl(i,j,k)   * frack
                    v_anl = vanl(i,j,k-1) * (1.-frack) 
     1                    + vanl(i,j,k)   * frack

                    u_rel_u = u_anl - ustorm(i,j)
                    v_rel_u = v_anl - vstorm(i,j)

                endif


!               Cross product of wind vectors on top and bottom of layer
                xprod = u_rel_l * v_rel_u - v_rel_l * u_rel_u       

!               Incremental area of hodograph
                area = .5 * xprod

!               Total area of hodograph
                area_sum = area_sum + area

                u_rel_l = u_rel_u
                v_rel_l = v_rel_u

            enddo ! k

            helicity(i,j) = area_sum * (-2.) 

            if(i .eq. 1 .and. j .eq. 1)then
                write(6,101)i,j,klow,khigh,helicity(i,j),ustorm(i,j)
     1                    ,(uanl(i,j,k),k=1,min(kmax,21))        
101             format(/4i4,f10.5,f8.3/7e11.3/7e11.3/7e11.3)
            endif

        enddo ! j
        enddo ! i

        return
        end
