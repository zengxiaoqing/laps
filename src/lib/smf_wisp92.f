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

      subroutine get_sfm_1d(nll,zcb,zctop,ityp,zlevel,
     1               plevel,tlevel,ql,qi,prob,mode)
      dimension zlevel(nll),plevel(nll),tlevel(nll),
     1          ql(nll),prob(nll),calw(200),cali(200),qi(nll)
cc   --------------------------------------------------------------
cc
cc         this is the streamlined version of the Smith-Feddes
cc         and Temperature Adjusted LWC calculation methodologies
cc         produced at Purdue University under sponsorship
cc         by the FAA Technical Center
cc
cc         currently, this subroutine will only use the Smith-
cc         Feddes and will only do so as if there are solely
cc         stratiform clouds present, however, it is very easy
cc         to switch so that only the Temperature Adjusted
cc         method is used
cc
cc         dilution by glaciation is also included, it is a
cc         linear function of in cloud temperature going from
cc         all liquid water at -10 C to all ice at -30 C
cc         as such the amount of ice is also calculated
cc
cc   --------------------------------------------------------------
      dz=100.0
      rv=461.5
      rair=287.04
      grav=9.81
      cp=1004.
      rlvo=2.5003e+6
      rlso=2.8339e+6
      dlvdt=-2.3693e+3
      eso=610.78
      c=0.01
      a1=8.4897
      b1=-13.2191
      c1=4.7295
      a2=10.357
      b2=-28.2416
      c2=8.8846
cc   -------------------------------------------------------------
cc
cc   calculate indices of cloud top and base
cc
cc   -------------------------------------------------------------
      nll_index = nll - 1
      do 100 ilevel=1,nll_index
        if(zlevel(ilevel).lt.zcb.and.
     1     zlevel(ilevel+1).gt.zcb) then
            kcb=ilevel
            kcb1=kcb+1
        endif
        if(zlevel(ilevel).lt.zctop.and.
     1     zlevel(ilevel+1).gt.zctop) then
            kctop=ilevel
            kctop1=kctop+1
        endif
100   continue
cc   -------------------------------------------------------------
cc
cc   obtain cloud base and top conditions
cc
cc   -------------------------------------------------------------
      delz=zlevel(kcb+1)-zlevel(kcb)
      delt=tlevel(kcb+1)-tlevel(kcb)
      cldbtm=delt*(zcb-zlevel(kcb))/delz+tlevel(kcb)
      tbar=(cldbtm+tlevel(kcb))/2.
      arg=-grav*(zcb-zlevel(kcb))/rair/tbar
      cldbp=plevel(kcb)*exp(arg)
      delz=zlevel(kctop+1)-zlevel(kctop)
      delt=tlevel(kctop+1)-tlevel(kctop)
      cldtpt=delt*(zctop-zlevel(kctop))/delz+tlevel(kctop)
cc   ---------------------------------------------------------------
cc
cc   calculate cloud lwc profile for cloud base/top pair
cc
cc   ---------------------------------------------------------------
        temp=cldbtm
        press=cldbp*100.0
        zbase=zcb
        nlevel=(zctop-zcb)/100.0
        if(nlevel.le.0) nlevel=1
        alw=0.0
        calw(1)=0.0
        nlm1=nlevel-1
        if(nlm1.lt.1) nlm1=1
        zht=zbase
        do 400 j=1,nlm1
          rl=rlvo+(273.15-temp)*dlvdt
          arg=rl*(temp-273.15)/273.15/temp/rv
          es=eso*exp(arg)
          qvs1=0.622*es/(press-es)
          rho1=press/(rair*temp)
          arg=-grav*dz/rair/temp
          p=press*exp(arg)
cc
cc   calculate saturated adiabatic lapse rate
cc
          des=es*rl/temp/temp/rv
          dtz=-grav*((1.0+0.621*es*rl/(press*rair*temp))/
     1        (cp+0.621*rl*des/press))
          zht=zht+dz
          press=p
          temp=temp+dtz*dz
          rl=rlvo+(273.15-temp)*dlvdt
          arg=rl*(temp-273.15)/273.15/temp/rv
          es2=eso*exp(arg)
          qvs2=0.622*es2/(press-es2)
          rho2=press/(rair*temp)
          alw=alw+(qvs1-qvs2)*(rho1+rho2)/2.
          calw(j+1)=alw
cc          write(6,9015) j,calw(j+1),zht
cc   ------------------------------------------------------------------
cc
cc   reduction of lwc by entrainment
cc
cc   ------------------------------------------------------------------
      ht=(zht-zbase)*.001
cc   ------------------------------------------------------------------
cc
cc                          skatskii's curve(convective)
cc
cc   ------------------------------------------------------------------
cc    if(ht.lt.0.3) then
cc      y=-1.667*(ht-0.6)
cc    elseif(ht.lt.1.0) then
cc      arg1=b1*b1-4.0*a1*(c1-ht)
cc      y=(-b1-sqrt(arg1))/(2.0*a1)
cc    elseif(ht.lt.2.9) then
cc      arg2=b2*b2-4.0*a2*(c2-ht)
cc      y=(-b2-sqrt(arg2))/(2.0*a2)
cc    else
cc      y=0.26
cc    endif
cc   ------------------------------------------------------------------
cc
cc                         warner's curve(stratiform)
cc
cc   ------------------------------------------------------------------
      if(ht.lt.0.032) then
        y=-11.0*ht+1.0
      elseif(ht.le.0.177) then
        y=-1.4*ht+0.6915
      elseif(ht.le.0.726) then
        y=-0.356*ht+0.505
      elseif(ht.le.1.5) then
        y=-0.0608*ht+0.2912
      else
        y=0.20
      endif
cc   -----------------------------------------------------------------
cc
cc                 calculate reduced lwc by entrainment
cc                             and dilution
cc
cc   -----------------------------------------------------------------
          ffrac=1.0
          if(temp.lt.263.15) then
            if(temp.gt.243.15) then
              ffrac=(temp-243.15)/20
            else
              ffrac=0.0
            endif
          endif
          if(ffrac.lt.0.0) ffrac=0.0
      tlwc=calw(j+1)
      calw(j+1)=calw(j+1)*y*1000.0*ffrac
      cali(j+1)=tlwc*y*1000.0-calw(j+1)
      if(cali(j+1).lt.0.0) cali(j+1)=0.0
 400  continue
cc
cc   -----------------------------------------------------------------
cc
cc   alternative calculation procedure using the observed or
cc   inferred in cloud temperature profile
cc
cc   ----------------------------------------------------------------
cc
      if(.true.) go to 455
        nlevel=(zctop-zcb)/100.0
        temp=cldbtm
        press=cldbp*100.0
        if(nlevel.le.0) nlevel=0
        alw=0.0
        calw(1)=0.0
        nlm1=nlevel-1
        if(nlm1.lt.1) nlm1=1
        dtdz=(cldtpt-cldbtm)/(zctop-zcb)
        zht=zbase
        do 450 j=1,nlm1
          rl=rlvo+(temp-273.15)*dlvdt
          arg=rl*(273.15-temp)/273.15/temp/rv
          es=eso*exp(arg)
          qvs1=0.622*es/(press-es)
          rho1=press/(rair*temp)
          arg=-grav*dz/rair/temp
          p=press*exp(arg)
          des=es*rl/temp/temp/rv
          dtz=-grav*((1.0+0.621*es*rl/(press*rair*temp))/
     1        (cp+0.621*rl*des/press))
          if(dtdz.lt.dtz) then
            dttdz=dtz-(dtdz-dtz)
          else
            dttdz=dtdz
          endif
          zht=zht+dz
          press=p
          temp=temp+dttdz*dz
          rl=rlvo+(273.15-temp)*dlvdt
          arg=rl*(temp-273.15)/273.15/temp/rv
          es2=eso*exp(arg)
          qvs2=0.622*es2/(press-es2)
          rho2=press/(rair*temp)
          alw=alw+(qvs1-qvs2)*(rho1+rho2)/2.
          if(alw.lt.0.0) alw=0.0
cc   ----------------------------------------------------------------
cc
cc   application of a simple linear glaciation
cc         all liquid T > -15 C
cc         partially liquid -15 C > T > -25 C
cc         all ice    T < -25 C
cc
cc   ---------------------------------------------------------------
          ffrac=1.0
          if(cldtpt.lt.258.15) then
            if(cldtpt.gt.248.15) then
              ffrac=(cldtpt-248.15)/10
            else
              ffrac=0.0
            endif
          endif
          if(ffrac.lt.0.0) ffrac=0.0
          calw(j+1)=alw*ffrac*1000.0
            write(6,9015) j,calw(j+1),zht
 450   continue
 455   continue
cc   -----------------------------------------------------------------
cc
cc                 obtain LAPS profile of LWCs
cc
cc   -----------------------------------------------------------------
      do 500 ip=1,nll
        if(zlevel(ip).lt.zcb.or.
     1     zlevel(ip).gt.zctop) then
            ql(ip)=0.0
            qi(ip)=0.0
        else
          ql(1)=calw(1)
          qi(1)=cali(1)
          do 600 j=2,nlevel
            zcloud=zcb+(j-1)*dz
            if(zcloud.ge.zlevel(ip)) then
              ql(ip)=(zlevel(ip)-zcloud+100.)*(calw(j)-
     1                    calw(j-1))/100.+calw(j-1)
              qi(ip)=(zlevel(ip)-zcloud+100.)*(cali(j)-
     1                    cali(j-1))/100.+cali(j-1)
              go to 650
            endif
 600      continue
 650      continue
        endif
 500  continue
cc   -----------------------------------------------------------------
cc
cc               write out file of lwc comparisons
cc
cc   -----------------------------------------------------------------
 9001 format(i7)
 9002 format(1x,3i2,1x,14f8.2,i2,i3)
 9004 format(1x,2e15.8,'ihr, imin, isec=',3(i8,1x))
 9005 format(2x,'Predicted LWC',8x,'Observed LWC',/)
 9014 format(1x,8e15.8)
 9015 format(1x,'j=',i3,'adiabatic lwc =',f12.6,'altitude =',f8.1)
      return
      end

