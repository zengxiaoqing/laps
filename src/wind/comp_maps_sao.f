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
        subroutine comp_maps_sao(u_bkg,v_bkg,ni,nj,nk
     1  ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sao,r_missing_data
     1  ,rms)

        real*4 grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
        real*4 grid_laps_wt(ni,nj,nk)
        dimension u_bkg(ni,nj,nk),v_bkg(ni,nj,nk) ! Model Background Data

        nobs = 0
        residualu = 0.
        residualv = 0.

        write(6,2)
2       format(
     1  /'Comparing Sao Velocities (passing QC) to Model Background'/
     1  1x,'   i   j   k      SAO                Model Background '
     1     ,'diff')

        do jl = 1,nj

          do il = 1,ni

!           if(grid_laps_wt(il,jl,nk/2) .eq. weight_sao)then ! for speed

              do k = 1,nk

                if(grid_laps_wt(il,jl,k) .eq. weight_sao
     1          .and. u_bkg(il,jl,k) .ne. r_missing_data )then
                  nobs = nobs + 1

                  diffu = u_bkg(il,jl,k) - grid_laps_u(il,jl,k)
                  diffv = v_bkg(il,jl,k) - grid_laps_v(il,jl,k)
                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2
                  write(6,101)il,jl,k
     1          ,grid_laps_u(il,jl,k),grid_laps_v(il,jl,k)
     1          ,u_bkg(il,jl,k),v_bkg(il,jl,k)
     1          ,-diffu,-diffv
101               format(1x,3i4,3(2x,2f7.1))

                endif

              enddo ! k

!           endif ! There is a SAO in this column

          enddo ! j
        enddo ! i


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
102     format(' RMS btwn Model Bkg & SAO (n,rmsu,rmsv,rms) = ',
     1     i4,3f5.1)

        return

        end

