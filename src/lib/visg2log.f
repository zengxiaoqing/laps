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
        subroutine visg2log(vis,ni,nj,smsng)
c
c       Routine to convert visibility from miles to log( miles ) for
c       the ni x nj laps grid.  If the visibility is zero, then set
c       the log( vis ) to -10 for the analysis.
c
c       01-14-92        P. Stamus
c
        real*4 vis(ni,nj)
c
        do j=1,nj
        do i=1,ni
          if(vis(i,j) .lt. 0.) then
            vis(i,j) = smsng
          elseif(vis(i,j) .eq. 0.) then
            vis(i,j) = -10.
          else
            vis(i,j) = log10( vis(i,j) )
          endif
        enddo !i
        enddo !j
c
        return
        end
