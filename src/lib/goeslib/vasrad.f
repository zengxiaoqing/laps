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
      real function vasrad(tau,temp,tsfc,kchan,nl,emiss)


c       routine received 5 aug 1993 (from wisconsin) as a replacement for
c       an earlier version that was out of date. and sent by mistake

c $ function vasrad(tau,temp,tsfc,kchan,nl,emiss)   (btr)
c $ vasrad - get vas radiance from temperature profile
c $ input:
c $     tau = (r) array of vas transmittances
c $     temp = (r) array of atmospheric temperatures
c $     tsfc = (r) temperature of surface
c $     kchan = (i) channel number
c $     nl = (i) level number of surface
c $     emiss = (r) surface emissivity
c $ output description:
c $     vas radiance
c $$ vasrad = sounder, vas, convert
      common/gde/gv(12),dv(12),ev(12)
      dimension tau(*),temp(*)
      t1=temp(1)
      b1=vplanc(t1,kchan)
      tau1=tau(1)
      rad=0.
      do 110 i=2,nl
      t2=temp(i)
      b2=vplanc(t2,kchan)
      tau2=tau(i)
      rad=rad+.5*(b1+b2)*(tau1-tau2)
      b1=b2
  110 tau1=tau2
      bs=vplanc(tsfc,kchan)
      rad=rad+emiss*bs*tau(nl)
      rad=rad+dv(kchan)
      rad=amax1(rad,.001)
      vasrad=rad
      return
      end
