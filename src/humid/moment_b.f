cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine moment_b(data,n,ave,adev,sdev,var,skew,curt,istatus)

c alternative moment program adapted from numerical recipies not to comflict
c with the other routine modifed and installed by Jim Edwards (moment.f)
c found under the satellite library area.


      integer n
      real adev,ave,curt,sdev,skew,var,data(n)
      integer j, istatus
      real p,s,ep
      istatus = 0 ! bad
      if(n.le.1) then
         write (6,*) 'must be at least 2 in moment comp'
         return
      endif
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      if(var.lt.0.) then
         write(6,*) 'negative VAR detected, taking absolute value'
         var = abs(var)
      endif
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
        skew = 0.0
        curt = 0.0
      write (6,*) 'Skew and Curt returned as zero, variance is zero'
        return
      endif
      istatus = 1 !good
      return
      end
