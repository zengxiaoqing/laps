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
       subroutine get_rirj_pol(latpt_in,
     &                         lonpt_in,
     &                         ny1,
     &                         lov,dx,dy,la1,lo1,
     &                         rj,
     &                         ri,
     &                         istatus)
c
c Program computes line and elem for satillite data available in line/elem
c from the ISPAN satellite data aquired in netCDF files (see subroutine
c readispan in calling program).
c
      implicit none

      real latpt_in
      real lonpt_in
      real u
      real v
      real ri
      real rj
      real pi 
      real du,dv
      real u_orig,v_orig
      integer ny1

      REAL*4      La1                
      REAL*4      Lo1                       
      REAL*4      Dx  
      REAL*4      Dy    
      REAL*4      LoV  

      integer*4 istatus
c
c ----------------------------- START --------------------------------------
c
c
      pi = acos (-1.0)
	
c get delta u and delta v in the polar grid

      call getdudv_pol(lov,dx/1000.,dy/1000.,la1,lo1,du,dv,
     &                 u_orig,v_orig)

c first get uv in polar grid

      call getuv_pol (pi/2.,lov*pi/180.,
     &                latpt_in*pi/180.,lonpt_in*pi/180.,
     &                u,v )

c compute ri, rj for polar stereographic data  

      call uv_ij (ny1,u_orig,v_orig,du,dv,u,v,ri,rj)
c
c note: ri and rj are returned!
c
      return
      end
