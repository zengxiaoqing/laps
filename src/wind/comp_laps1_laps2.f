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
        subroutine comp_laps1_laps2(uf,vf,u,v,ni,nj,nk,rms)

        include 'windparms.inc' 

        dimension u(ni,nj,nk),v(ni,nj,nk)
        dimension uf(ni,nj,nk),vf(ni,nj,nk)

        call get_r_missing_data(r_missing_data, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in comp_laps1_laps2'
            stop
        endif

        nobs = 0
        residual = 0.
        nprint = 0

        write(6,2)
2       format(/'      Comparing Both LAPS analyses'
     1       /'   i   j   k     u1     v1         u2     v2  '
     1       ,'      du    dv')

        do k = 1,nk
        do j = 1,nj
        do i = 1,ni

            if(        u(i,j,k) .ne. r_missing_data
     1          .and. uf(i,j,k) .ne. r_missing_data)then
                nobs = nobs + 1

                diffsq = (u(i,j,k) - uf(i,j,k))**2 +
     1                   (v(i,j,k) - vf(i,j,k))**2
                residual = residual + diffsq
                if((nobs/800) * 800 .eq. nobs .or. diffsq .gt. 20.)then
                  if(nprint .lt. 200)then
                    write(6,11)i,j,k,uf(i,j,k),vf(i,j,k)
     1                        ,u(i,j,k),v(i,j,k)
     1                        ,u(i,j,k) - uf(i,j,k)
     1                        ,v(i,j,k) - vf(i,j,k)
11                  format(1x,3i4,2f7.1,3x,2f7.1,3x,2f7.1)
                    nprint = nprint + 1
                  endif
                endif

            endif

        enddo ! i
        enddo ! j
        enddo ! k

        if(nobs .gt. 0)then
            rms = sqrt(residual/nobs)
        else
            rms = 0.
        endif

        write(6,*)' RMS between First and Second & LAPS anal = ',nobs,rm
     1s
        write(15,*)' RMS between First and Second & LAPS anal = ',nobs,r
     1ms

        return

        end

