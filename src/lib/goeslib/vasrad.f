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
