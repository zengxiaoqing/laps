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
      subroutine draw_county_map(sw,ne,jproj,polat,polon,rrot,jdot
     1                                             ,icol_sta,icol_cou)
c
c
      real*4 sw(2),ne(2),pl1(2),pl2(2),pl3(2),pl4(2),
     +       polat,polon,rrot

c
      integer*4 jproj,jjlts,jgrid,jus,jdot,ier
c
      COMMON/SUPMP9/DS,DI,DSRDI

!     DI = 50.
!     polat=90.

      rrot=0.
      pl1(1)=sw(1)
      pl2(1)=sw(2)
      pl3(1)=ne(1)
      pl4(1)=ne(2)
      jjlts=-2
      jgrid=0


!     Plot Counties
      jus=-4

      if(jdot .eq. 1)then
          call gsln(3) ! Dotted
      else
          call gsln(1) ! Solid
      endif

      call setusv_dum(2HIN,icol_cou)
      call supmap(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,jjlts,
     +            jgrid,jus,jdot,ier)
      call sflush

!     Plot States From Counties
      jus=-8
      call gsln(1)
      call setusv_dum(2HIN,icol_sta)
      call supmap(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,jjlts,
     +            jgrid,jus,jdot,ier)

!     Plot Continents
      jus=-1
      call gsln(1)
      call setusv_dum(2HIN,icol_sta)
      call supmap(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,jjlts,
     +            jgrid,jus,jdot,ier)

      return
      end
