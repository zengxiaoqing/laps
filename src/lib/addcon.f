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
c
c
        subroutine addcon(a,con,result,imax,jmax)
c
cdoc    Routine to add the constant 'con' to array 'a' and put the
cdoc    result into array 'result'.
c
        real a(imax,jmax), result(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          result(i,j) = a(i,j) + con
        enddo !i
        enddo !j
c
        return
        end


        subroutine addcon_miss(a,con,result,imax,jmax)
c
cdoc    Routine to add the constant 'con' to array 'a' and put the
cdoc    result into array 'result'. This takes account of 'r_missing_data'.
c
        real a(imax,jmax), result(imax,jmax)
c
        call get_r_missing_data(r_missing_data,istatus)

        do j=1,jmax
        do i=1,imax
          if(a(i,j) .ne. r_missing_data .and. 
     1       con    .ne. r_missing_data       )then
              result(i,j) = a(i,j) + con
          else
              result(i,j) = r_missing_data
          endif
        enddo !i
        enddo !j
c
        return
        end
