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
        subroutine comp_laps_sao(u,v,ni,nj,nk,n_sao_obs
     1          ,sao_i,sao_j,sao_k,sao_u
     1          ,sao_v,N_SAO,r_missing_data,rms)

        integer sao_i(N_SAO) ! X sao coordinates
        integer sao_j(N_SAO) ! Y sao coordinates
        integer sao_k(N_SAO) ! Z sao coordinates
        real    sao_u(N_SAO) ! u sao component
        real    sao_v(N_SAO) ! v sao component

        dimension u(ni,nj,nk),v(ni,nj,nk)

          nobs = 0
          residualu = 0.
          residualv = 0.

          write(6,2)
2         format(/'        Comparing LAPS to Saos,   '/
     1  1x,'   n   i   j   k       sao        ',
     1          '  LAPS          diff ')

          do n = 1,n_sao_obs

            i = sao_i(n)
            j = sao_j(n)
            k = sao_k(n)
            u_sao = sao_u(n)
            v_sao = sao_v(n)

            if(k .le. nk .and. i .gt. 0 .and. j .gt. 0 .and. i .le. ni .
     1and.
     1       j .le. nj)then

              if(u(i,j,k) .ne. r_missing_data
     1   .and. v(i,j,k) .ne. r_missing_data)then
                nobs = nobs + 1
                diffu = u_sao - u(i,j,k)
                diffv = v_sao - v(i,j,k)
                residualu = residualu + diffu ** 2
                residualv = residualv + diffv ** 2
                write(6,101)n,i,j,k,u_sao,v_sao,
     1                  u(i,j,k),v(i,j,k),diffu,diffv
101                 format(1x,4i4,3(1x,2f6.1))

              else
                write(6,*)' NO LAPS WINDS AT THIS SAO LOCATION'

              endif

            endif ! k

          enddo ! n

          if(nobs .gt. 0)then
              rmsu = sqrt(residualu/nobs)
              rmsv = sqrt(residualv/nobs)
              rms  = sqrt(rmsu**2 + rmsv**2)
          else
              rmsu = 0.
              rmsv = 0.
              rms = 0.
          endif

          write(6 ,102)nobs,rmsu,rmsv,rms
          write(15 ,102)nobs,rmsu,rmsv,rms
102       format(' RMS between LAPS & sao     (n,rmsu,rmsv,rms) = ',
     1     i4,3f5.1)


        return

        end

