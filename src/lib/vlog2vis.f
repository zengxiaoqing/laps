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
        subroutine vlog2vis(lvis,vis,ni,nj)
c
c       Routine to convert gridded visibility from log( miles ) to
c       miles.  If the analyzed log (vis) is less that -2, set the
c       visibility to zero.
c
c       01-14-92        P. Stamus
c
        real*4 vis(ni,nj), lvis(ni,nj)
c
        do j=1,nj
        do i=1,ni
          if(lvis(i,j) .le. -1.5) then
            vis(i,j) = 0.
          else
            vis(i,j) = 10. ** ( lvis(i,j) )
          endif
        enddo !i
        enddo !j
c
        return
        end
