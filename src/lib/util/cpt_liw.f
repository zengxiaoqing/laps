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

        subroutine cpt_liw(lifted,w_2d,imax,jmax,liw)

!       Note, these values are not logarithmically scaled
!       Li * Omega is returned in Pa K/s

        real*4 w_2d(imax,jmax) ! Omega (Pa/s)
        real*4 lifted(imax,jmax) ! Lifted Index (K)
        real*4 liw(imax,jmax) ! Pa K/s

        do j = 1,jmax,1
        do i = 1,imax,1

            if(lifted(i,j) .lt. 0. .and. w_2d(i,j) .lt. 0.)then
                liw(i,j) = lifted(i,j) * w_2d(i,j)

            else
                liw(i,j) = -1e-6

            endif

        enddo ! j
        enddo ! i

        return
        end

