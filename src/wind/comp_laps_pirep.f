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
        subroutine comp_laps_pirep(u,v,ni,nj,nk,r_missing_data
     1  ,weight_pirep,grid_laps_u,grid_laps_v,grid_laps_wt,rms)

        dimension u(ni,nj,nk),v(ni,nj,nk)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)
        real grid_laps_wt(ni,nj,nk)

        nobs = 0
        residualu = 0.
        residualv = 0.

        write(6,2)
2       format(/'      Comparing Pirep Velocities to LAPS'/
     1  6x,'   i   j   k    pirep              laps             diff ')       

        do k = 1,nk
        do j = 1,nj
        do i = 1,ni

            if(       grid_laps_wt(i,j,k) .eq. weight_pirep
     1          .and. u(i,j,k)            .ne. r_missing_data )then
                nobs = nobs + 1

                diffu = u(i,j,k) - grid_laps_u(i,j,k)
                diffv = v(i,j,k) - grid_laps_v(i,j,k)
                residualu = residualu + diffu ** 2
                residualv = residualv + diffv ** 2

                if(nobs .le. 200 .OR. nobs .eq. (nobs/10)*10)then
                    write(6,101)nobs,i,j,k
     1              ,grid_laps_u(i,j,k),grid_laps_v(i,j,k),
     1              u(i,j,k),v(i,j,k),diffu,diffv
101                 format(1x,i5,3i4,3(2x,2f7.1))
                endif

            endif

        enddo ! i
        enddo ! j
        enddo ! k


        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            rmsv = sqrt(residualv/nobs)
        else
            rmsu = 0.
            rmsv = 0.
        endif

        rms  = sqrt(rmsu**2 + rmsv**2)

        write(6 ,102)nobs,rmsu,rmsv,rms
        write(15 ,102)nobs,rmsu,rmsv,rms
102     format(' RMS between LAPS & Pireps   (n,rmsu,rmsv,rms) = ',
     1     i4,3f5.1)

        return

        end

