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

