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
        subroutine move(a,b,imax,jmax)
c
c.....  Routine to move array 'a' into array 'b'.
c
        real a(imax,jmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end
c
c-----------------------------------------------------
c
        subroutine move_i(a,b,imax,jmax)
c
c.....  Routine to move array 'a' into array 'b'.
c
        integer a(imax,jmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end

c
c-----------------------------------------------------
c
        subroutine move_3d(a,b,imax,jmax,kmax)
c
c.....  Routine to move (copy) array 'a' into array 'b'.
c
        real a(imax,jmax,kmax), b(imax,jmax,kmax)
c
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
          b(i,j,k) = a(i,j,k)
        enddo !i
        enddo !j
        enddo !k
c
        return
        end
c
c
        subroutine move_2dto3d(a,b,index,imax,jmax,kmax)
c
c.....  Routine to move (copy) the 2d array 'a' into one level
c.....  of the 3d array 'b'.  The level is defined by 'index'.
c
c       Original:  P. Stamus  NOAA/FSL  15 Apr 1997
c
        real a(imax,jmax), b(imax,jmax,kmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j,index) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end
c
c
        subroutine move_3dto2d(a,index,b,imax,jmax,kmax)
c
c.....  Routine to move (copy) one level of the 3d array 'a' into 
c.....  the 2d array 'b'.  The level is defined by 'index'.
c
c       Original:  P. Stamus  NOAA/FSL  15 Apr 1997
c
        real a(imax,jmax,kmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j,index)
        enddo !i
        enddo !j
c
        return
        end
