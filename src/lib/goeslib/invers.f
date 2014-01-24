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

        subroutine filter(vdat,nspin,erm,msam,tsfs,ldetr)
c
c ****  tsfs is surface air temp. (from surface analysis)
c
c      dimension vdat(13),nspin(13),erm(13),lbuf(33)
      dimension vdat(13),nspin(13),erm(13)
c
c ****  added icloud to filt common block.
      common/filt/tdat(13,11,11),iflag(2,11,11),icloud
      common/dbug/kbug
      common/gde/gv(12),dv(12),ev(12)
c
c ****  added nmove to size common block.
      common/size/nbxs,nmove
      data vmisg/999999./
      minsam=nbxs*nbxs/20       !changed
      minsam=min0(minsam,2)     !line moved
      nbxm=minsam+1
      minsam=max0(minsam,1)
      vmax=1000.                !new
      vmin=tsfs-25.             !new
      itr=0                     !new
50    msam=0                    !added line number
      itr=itr+1                 !new
      vsfs=tsfs                 !changed tsfs-10. --> tsfs
      v=vmin                    !changed 0. --> vmin
      vmin5=vmisg
      vtdif=0.                  !new
      do 100 j=1,nbxs
      do 100 i=1,nbxs
      v5=tdat(5,i,j)
      if(v5.eq.vmisg) go to 100
      v7=tdat(7,i,j)            !new
      if(v7.eq.vmisg) go to 100 !new
      v8=tdat(8,i,j)
      if(v8.eq.vmisg) go to 100
      if(v8.ge.vmax) go to 100  !new
      if(v8.lt.v) go to 100
      v=v8
      vmin7=v7                  !new
      vmin5=v5
      if(tdat(12,i,j).eq.vmisg) go to 100       !new
      vtdif=tdat(12,i,j)-tdat(8,i,j)            !new
  100 continue
c
c-----deleted next section.-----------------------------------------------------
c     if(v.lt.vsfs) return
c     if(vmin5.eq.vmisg) return
c     tmin=v-2.
c     vmin5=vmin5-1.
c
c-----added new section.--------------------------------------------------------
c
      vmax=v
      if(vmin5.eq.vmisg) return
c
c ****  set min temp to 3 degrees less than warmest channel 8.
c
      tmin=v-3.
      vmin5=vmin5-1.
c
c ****  set up reflected sunlight test using vmintd.
c
      vmintd=vtdif+2.
      vmintd=amin1(5.,vmintd)
      icloud=0
c
c ****  consider channel 12 emissivity.
c ****  observation suggests that inversions will never give a -4..maybe
c ****  note that we repeat this test after filtering.
c
      if(vtdif.lt.-4.) go to 103
c
c ****  avoid reflected sunlight in check below.
c
      if(vtdif.gt.0) vtdif=0
c
c ****  make split window adjustment to estimate 11 micron from tsfc.
c
      spw=v-vmin7
      if (spw.lt.0)spw=0.
      vsfs=vsfs-spw
      vtmin=3.
      if(v.gt.270.) vtmin=5.
c
c ****  check warmest window against adjusted surface.
c
      vchk=vsfs-v+vtdif
c
c ****  note that an inversion will allow greater tsfc-v8 discrepancy.
c
      if(vck.lt.vtmin) go to 105
103   tmin=tmin-2.
c
c ****  reduce limit on cloudy.
c ****  check rediculously cold.
c
      if(tmin.lt.vmin) return
      vmin5=vmin5-1.
      icloud=1
105   continue
      tmin=amin1(tmin,tsfs)
c
c-----end new section.----------------------------------------------------------
c
c               02-may-84       ibm to vax fortran change
c     if(kbug.ne.0)call enkode('(132x,t1,"begin filter, tmin =",f7.2/)',
c    * lbuf,tmin)
        if (kbug .ne. 0) write (6,5000) tmin
5000    format (1x, 'begin filter, tmin =', f7.2)
c * produce large-detector channel 8
      do 160 j=1,nbxs
      ll=j-1
      ll=max0(ll,1)
      lh=j+1
      lh=min0(lh,nbxs)
      do 150 i=1,nbxs
      iflag(1,i,j)=0
      iflag(2,i,j)=0
      tdat(13,i,j)=vmisg
c     if(tdat(6,i,j).eq.vmisg) go to 150        !line deleted
      if(tdat(8,i,j).eq.vmisg) go to 150
      iel=i-1
      iel=max0(iel,1)
      ieh=i+1
      ieh=min0(ieh,nbxs)
      v=vmisg
      nv=0
      do 120 il=ll,lh
      do 110 ie=iel,ieh
      tdat6=tdat(6,ie,il)
      if(tdat6.eq.vmisg) go to 110
      if(tdat6.lt.tdat(6,i,j)) go to 110
      tdat8=tdat(8,ie,il)
      if(tdat8.eq.vmisg) go to 110
      if(tdat8.lt.tmin) go to 110
      v=amin1(v,tdat8)
      nv=nv+1
  110 continue
  120 continue
      if(nv.lt.4) go to 140
      tdat(13,i,j)=v
      if(kbug.eq.0) go to 140
c
c               02-may-84       ibm to vax fortran change
c     call enkode('(" j, i, l-d tbb(8) =",2i3,f7.2/)',lbuf,j,i,v)
        write (6, 5010) j, i, v
5010    format (1x, ' j, i, l-d tbb(8) =', 2i3, f7.2)
  140 continue
c * cloud filter
      if(tdat(8,i,j).eq.vmisg) go to 150        !removed .or. .lt. tmin
      df=0.                                     !new
      if(tdat(12,i,j).ne.vmisg) df=tdat(12,i,j)-tdat(8,i,j)     !new
      if(tdat(8,i,j).lt.tmin) go to 150         !new
      if(tdat(8,i,j).gt.vmax) go to 150         !new
c
c ****  accept long wave.
c
      iflag(1,i,j)=1
      if (df.gt.vmintd)go to 150                !new
      if(tdat(13,i,j).eq.vmisg) go to 150
      iflag(2,i,j)=1
  150 continue
  160 continue
c
      do 250 k=1,12
      num=0
      sum=0.
      vdat(k)=vmisg
      erm(k)=vmisg
      do 240 j=1,nbxs
      do 230 i=1,nbxs
      tdk=tdat(k,i,j)
      if(tdk.eq.vmisg) go to 230
      if(k.lt.3) go to 220
c     if(k.eq.6.or.k.eq.12) go to 210           !line deleted
c
c ****  lump channel 6 in with the longwave.
c
      if(k.eq.12)go to 210                      !new
      if(k.gt.4.and.k.lt.9) go to 200
c
c ****  filter water vapor channels.
c
180   td5=tdat(5,i,j)                           !added line number
      if(td5.eq.vmisg) go to 230
      if(td5.gt.vmin5) go to 220
  200 if(iflag(1,i,j).eq.0) go to 230
      go to 220
  210 if(iflag(2,i,j).eq.0) go to 230
  220 num=num+1
      sum=sum+tdk
  230 continue
  240 continue
      if(sum.eq.0.) go to 250
      if(num.eq.0) go to 250
      if(nspin(k).eq.0) go to 250
      sum=sum/float(num)
      if(sum.lt.180..or.sum.gt.330.) go to 250
      vdat(k)=sum
      dbdtbb=vdbdtb(vdat(k),k)
      spfac=1.
      if(ldetr.ne.0) go to 245
      if(k.lt.3.or.k.eq.6.or.k.gt.8) spfac=0.25
  245 spins=spfac*num*nspin(k)
      erm(k)=(ev(k)/dbdtbb)/sqrt(spins)
      if(k.ne.8) go to 250
c     if(num.lt.minsam) return          !deleted line
      msam=num
  250 continue
      if(erm(7).ne.vmisg) erm(7)=amax1(erm(7),0.5)
      if(erm(8).ne.vmisg) erm(8)=amax1(erm(8),0.5)
c
c-----new section.--------------------------------------------------------------
c
c ****  if we didn't get a channel 12 sample we have been looking at
c         clouds; either thin cirrus (at night) or reflected solar.
c
      if(vdat(12).eq.vmisg) go to 253
c
c ****  repeat cloud tests.
c
      icloud=0
c
c ****  consider channel 12 emissivity.
c
      vtdif=vdat(12)-vdat(8)
      if(vtdif.lt.-4.) go to 252
      if(vtdif.gt.0) vtdif=0
      spw=vdat(8)-vdat(7)
      if (spw.lt.0.)spw=0.
      vsfs=tsfs-spw
      vtmin=3.
      if(v.gt.270.) vtmin=5.
      vck=vsfs-vdat(8)+vtdif
      if(vck.lt.vtmin) go to 255
252   tmin=tmin-2.
      if(tmin.lt.vmin) msam=0
253   icloud=1
c
c-----end new section.----------------------------------------------------------
c
255   if(kbug.eq.0) return              !added line number
c
c               02-may-84       ibm to vax fortran change
c     call enkode('(" filter complete, chan-8 sample =",i4/)',lbuf,msam)
        write (6, 5020) msam
5020    format (1x, ' filter complete, chan-8 sample =', i4)
      return
      end
      subroutine qvtwr(tbo,ero,toto,tskin,tsta,wsta,lsta,nsat,ifail)
      common/atmos/pre(40),tem(40),wmx(40)
      common/dbug/kbug
      common/func/phi(40,9),xit(9,12),ers(12),tbc(10),nch,nft,nfw
c
c ****  add ipass to gam common block
      common/gam/gamma,est,esw,ipass
      common/nav/vlat,vlon,vzen,szen,il,ie,iras,ipic,itime,jtime,jday
c
c ****  new common blocks.
      common/cloud/ms,tcld,wcld
      common/radtra/tau,dbdt
      common/filt/tdat(13,11,11),iflag(2,11,11),icloud
      real*8 xtx(9,9),xiv(9,9)
c      dimension tbo(*),ero(*),dtb(12),lbuf(33)
      dimension tbo(*),ero(*),dtb(12)
      dimension alpha(9,12),phs(40,3),phis(40,3),coef(9)
      dimension tau(40),taus(40,10),dbdt(40),pwv(40),delp(40),dlnp(40)
c
c ****  change item(4) --> item(5) & iwat(2) --> iwat(3)
      dimension tauh(40),dtdu(40),pwvh(40),wmxh(40),item(5),iwat(3)
      dimension ts(40),ws(40)                   !new
      data nt/5/,nw/3/
      data item/2,4,5,7,3/                      !added band 3
      data iwat/7,9,10/                         !added band 10
      data init/1/,nchan/10/,nl/40/,ndim/9/     !added ndim
      data grav/980.665/,vmisg/999./


c       vax to ibm to cray fortran change, initialize var "bc".  db/4/13/94
        bc = 0.0
c

c
c
c               02-may-84       ibm to vax fortran change
c     if(kbug.ne.0) call enkode('("entering qvtwr ..."/)',lbuf)
        if (kbug .ne. 0) write (6, 5000)
5000    format (1x, 'entering qvtwr ...')
      if(init.eq.0) go to 120
      nch=nchan
c
c ****  temporary band-aids.   07-mar-86
c ****  removed.               28-apr-86
c
c       ms=lsta
c       icloud=0
c
        if (kbug .ne. 0) then
                write (6, 7000) icloud
7000            format(1x,' icloud =',i2)
        endif
c
c ****  set up indices.
c
      nft=nt
      nfw=nw
      ntt=nt-1
      nwt=nw-1
      nwp1=nw+1
      ntot=nw+nt+1
      ntm1=ntot-1
      delp(1)=0.
      do 100 i=2,nl
      dlnp(i)=alog(pre(i)/pre(i-1))
  100 delp(i)=pre(i)-pre(i-1)
      init=0
c
c
c            1 5-aug-84 debug
c
c 120 ifail=1
120     ifail = 0
      do 140 k=1,nchan
      if(tbo(k).ge.vmisg) ero(k)=1.e9
  140 continue
c
c ****  new
c ****  on initial entry, modify guess mixing ratio from sfc. to 100mb.
c
      if (ipass .eq. 1) then
                ls=lsta
                cor=wsta/wmx(ls)
                do 160 i=21,ls
                        wmx(i)=wmx(i)*cor
  160           continue
      endif
c
c ****  get first estimate of transmittance for all channels.
c
      do 180 k=1,nchan
                kc=k
                call vastau(tem,wmx,toto,vzen,taus(1,k),nsat,kc)
  180 continue
c
c ****  form normalized water vapor basis functions from 'iwat' channels.
c
      do 240 j=1,nw                             !changed nwt --> nw
                k=iwat(j)
                do 200 i=1,nl
                        tau(i)=taus(i,k)
  200           continue
c
c ************  uses phs array as scratch before weighting with
c                 mixing ratio and integrating.
c
                phs(1,j)=0.
                psmax=0.
                do 220 i=2,nl
                        dt=tau(i-1)-tau(i)
                        phs(i,j)=dt/dlnp(i)
c
c ********************  new
c ********************  set function to zero below cloud level.
c
                        if (i .gt. ms) phs(i,j)=0.
                        psmax=amax1(psmax,phs(i,j))
  220           continue
                do 230 i=2,nl
                        phs(i,j)=phs(i,j)/psmax
  230           continue
  240 continue
c
c ****  form normalized temperature weighting functions.
c
      do 300 j=1,nt                             !changed ntt --> nt
                k=item(j)
                do 260 i=1,nl
                        tau(i)=taus(i,k)
  260           continue
c
c ************  index temperature basis functions after water vapor.
c
                jj=j+nw
                phi(1,jj)=0.
                psmax=0.
                do 280 i=2,nl
                        dt=tau(i-1)-tau(i)
                        phi(i,jj)=dt/dlnp(i)
                        if (i .gt. ms) phi(i,jj)=0.     !new
                        psmax=amax1(psmax,phi(i,jj))
  280           continue
                do 290 i=2,nl
                        phi(i,jj)=phi(i,jj)/psmax
  290           continue
  300 continue
c
c ****  new
c ****  for first pass do not use channels which conflict (3 and 10).
c
      if (ipass .eq. 1) then
                do 320 i=1,nl
                        phs(i,nw)=1.0
                        phi(i,ntm1)=1.0
  320           continue
      endif
c
c ****  complete integration of water band function and store in phis.
c
      do 380 j=1,nw
                do 340 i=1,nl
                        ph=phs(i,j)
                        phi(i,j)=ph*wmx(i)/1000.
  340           continue
                phis(1,j)=0.
                sum=0.
                ph1=phi(1,j)
                do 360 i=2,nl
                        ph2=phi(i,j)
                        sum=sum+0.5*(ph1+ph2)*delp(i)
                        phis(i,j)=sum
                        ph1=ph2
  360           continue
  380 continue
c
c ****  save the guess profiles.
c
      do 400 i=1,nl
                ts(i)=tem(i)                    !new
                ws(i)=wmx(i)                    !new
                wmxh(i)=0.5*wmx(i)
  400 continue
c
c ****  set up precip water profiles for observed and 0.5*observed.
c         in order to solve for dtau/du.
c
      call precw(pre,wmxh,pwvh,nl)
      call precw(pre,wmx,pwv,nl)
c
c               02-may-84       ibm to vax fortran change
c     if(kbug.ne.0) call enkode('("tskin =",f7.2/)',lbuf,tskin)
        if (kbug .ne. 0) write (6,5010) tskin
5010    format (1x, 'tskin =', f7.2)
c
c ****  new
c ****  set up reflected sunlight check.
c
c     t12chk=tbo(12)-tskin              !line deleted 25-apr-86
c
c ****  new
c ****  initialize entire alpha array.
c
      do 470 k=1,12
      do 470 j=1,9
                alpha(j,k)=0.
  470 continue
c
      do 560 k=1,nchan
                kc=k
c
c ************  fill tauh as transmittance for 1/2 precipitable water.
c
                if (k .ge. 3) then
                        call vastau(tem,wmxh,toto,vzen,tauh,nsat,kc)
                endif
                do 420 i=1,nl
                        tau(i)=taus(i,k)
  420           continue
                is=ls                                   !new
                if (ls .gt. ms) is=ms                   !new
c
c ************  calculate brightness t and partial derivatives needed for
c                 basis f.
c
                call vasrte(tau,tem,tskin,rad,dbdt,tbc(k),dbdtbb,kc,is)
c                                                       -changed ls to is
                dtb(k)=tbo(k)-tbc(k)
                adtb=abs(dtb(k))
c
c ************  put lower limit on error.
c
                erkc=amax1(ero(k),0.25)
                if(k.ne.6) go to 440
                ers6=ero(k)                             !new
                if (k .eq. 6) go to 440                 !added 22-aug-1986
                                                        !do not check ch6.
c
c ************  save channel 6 for 2nd pass reevaluation
c                 in event it is tossed on first pass.
c                 default option is to not check channel 6 consistency.
c
                if (tbo(k) .ge. vmisg) goto 440         !new
                if(adtb.lt.0.5) go to 440               !changed 1.0 --> 0.5
c
c ************  accept channel 6 if it agrees with channel 5.
c
                if (abs(dtb(6)-dtb(5)) .lt. 0.5) goto 440       !new
c
c ************  look at cloud flag from subroutine filter and reject
c                 channel 6 if on.
c
                if (icloud .eq. 1) goto 430                     !new
                if(adtb.lt.abs(dtb(5))) go to 440
430             adtb=1000.                                      !new
                erkc=adtb
  440           ers(k)=1./erkc
                ero(k)=erkc
c
c               02-may-84       ibm to vax fortran change
c     if(kbug.ne.0)
c    *call enkode('("chan ",i2,", tbo =",f7.2,", tbc =",f7.2,", dtb =",
c    * f6.2,", ero =",f6.2/)',lbuf,kc,tbo(k),tbc(k),dtb(k),ero(k))
        if (kbug .ne. 0) write (6, 5020) kc, tbo(k), tbc(k),
     1  dtb(k), ero(k)
5020    format (1x, 'chan ', i2, ', tbo =', f7.2, ', tbc =', f7.2,
     1        ', dtb =', f6.2, ', ero =', f7.2)
c
c ************  new
c ************  set up array of partial derivatives.
c
                do 450 i=1,ls
                        dbdt(i)=dbdt(i)/dbdtbb
  450           continue
                if (k .lt. 3) goto 510
c
c ************  calculate dtau/du
c
                do 460 i=1,ls
                        dtdu(i)=0.
                        dt=tau(i)-tauh(i)
                        if(dt.eq.0.) go to 460
                        du=pwv(i)-pwvh(i)
                        if (abs(du) .lt. 1.e-5) goto 460        !new
                        dtdu(i)=dt/du
  460           continue
c
c ************  complete basis function for water vapor.
c
                do 500 j=1,nw
                        sum=0.
                        phd1=phis(1,j)*dbdt(1)*dtdu(1)
                        tem1=tem(1)
                        do 480 i=2,is                   !changed ls --> is
                                phd2=phis(i,j)*dbdt(i)*dtdu(i)
                                tem2=tem(i)
                                sum=sum+.5*(phd1+phd2)*(tem2-tem1)
                                phd1=phd2
                                tem1=tem2
  480                   continue
                        alpha(j,k)=sum
  500           continue
c
  510           do 540 j=nwp1,ntm1
                        sum=0.
                        phd1=phi(1,j)*dbdt(1)
                        tau1=tau(1)
                        do 520 i=2,is                   !changed ls --> is
                                phd2=phi(i,j)*dbdt(i)
                                tau2=tau(i)
                                sum=sum+.5*(phd1+phd2)*(tau2-tau1)
                                phd1=phd2
                                tau1=tau2
  520                   continue
                        alpha(j,k)=sum
  540           continue
c
c ************  set 0th order basis function.
c
                alpha(ntot,k)=dbdt(is)*tau(is)          !changed ls --> is
  560 continue
c
c ****  add "observation" of surface mixing ratio and compute functions.
c
      k=nchan+1
      dtb(k)=wsta-wmx(ls)
      ers(k)=1./(esw*wsta)
      do 580 j=1,nw
                alpha(j,k)=grav*phi(ls,j)
  580 continue
      do 600 j=nwp1,ntot
                alpha(j,k)=0.
  600 continue
c
c ****  add observation of surface temperature and calculate
c         contribution for temperature functions.
c
      k=k+1
      dtb(k)=tem(ls)-tsta
      ers(k)=1./est
      do 620 j=nwp1,ntm1
                alpha(j,k)=phi(ls,j)
  620 continue
      do 640 j=1,nw
                alpha(j,k)=0.
  640 continue
      alpha(ntot,k)=0.
      neqn=k
c
c ****  new
c ****  if cloudy, replace two surface equations with cloud level terms.
c ****  note: on first pass these will have no effect (wcld=wmx(ms)).
c
      if (ms .lt. ls) then
                k=k-1
                dtb(k)=wcld-wmx(ms)
                ers(k)=1./(esw*wcld)
                do 643 j=1,nw
                        alpha(j,k)=grav*phi(ms,j)
643             continue
                do 646 j=nwp1,ntot
                        alpha(j,k)=0.
646             continue
                k=k+1
                dtb(k)=tem(ms)-tcld
                ers(k)=1./est
                do 648 j=nwp1,ntm1
                        alpha(j,k)=phi(ms,j)
648             continue
                do 650 j=1,nw
                        alpha(j,k)=0.
650             continue
                alpha(ntot,k)=0.
      endif
c
c ****  multiply both sides of equation by 1/error.
c
      do 680 k=1,neqn
                es=ers(k)
                do 660 j=1,ntot
                        alpha(j,k)=alpha(j,k)*es
  660           continue
                dtb(k)=dtb(k)*es
  680 continue
c
c               07-jun-84       debug
c
        if (kbug .ne. 0) write (6, 5050) ((alpha(iii, jjj), iii = 1, 9),
     1                                               jjj = 1, 12)
5050    format (/, 1x, 'alpha', /, (9(1x, f10.4)))
      call solvex(gamma,alpha,xit,xtx,xiv,ndim,neqn)    !changed ntot --> ndim
      do 690 j=1,ntot
                sum=0.
                do 685 k=1,neqn
                        sum=sum+xit(j,k)*dtb(k)
  685           continue
                coef(j)=sum
  690 continue
c
c               07-jun-84       debug
c
        if (kbug .ne. 0) write (6, 5060) (coef(iii), iii = 1, 9)
5060    format (/, 1x, 'coef', /, 9(1x, f10.4))
c
      do 720 i=21,is                                    !changed ls --> is
                sum=0.
                do 700 j=1,nw
                        sum=sum+coef(j)*phi(i,j)
  700           continue
                wmx(i)=wmx(i)+grav*sum
  720 continue
      do 760 i=1,is                                     !changed ls --> is
                sum=0.
                do 740 j=nwp1,ntm1
                        sum=sum+coef(j)*phi(i,j)
  740           continue
                test=abs(sum)
c
c            1 3-aug-83 set error value
c
        ifail = 1
c
                if(i.gt.26.and.test.gt.10.) return
                tem(i)=tem(i)-sum
 760  continue
      test=abs(coef(ntot))
c
c            1 3-aug-84 set error value
c
        ifail = 2
      if(test.gt.10.) return
      if (bc .eq. 0) tskin=tskin+coef(ntot)
c
c ****  new
c ****  below cloud, blend profile with surface data.
c
      if (is .ne. ls) then
                dp=pre(is)-pre(ls)
                dtdp=(tem(is)-ts(is))/dp
                dwdp=(wmx(is)-ws(is))/dp
                do 765 i=is,ls
                        dp=pre(i)-pre(ls)
                        tem(i)=ts(i)+dtdp*dp
                        wmx(i)=ws(i)+dwdp*dp
  765           continue
      endif
c
      do 780 i=ls,nl
                tem(i)=tem(ls)
                wmx(i)=wmx(ls)
  780 continue
      do 800 i=21,nl
                wmax=0.95*wsat(pre(i),tem(i))           !added 0.95*
                wmin=wmax*0.0002*pre(i)         !changed 0.05 --> 0.0002*pre(i)
                wlim=2.0*wmax                   !changed 1.5 --> 2.0
c
c            1 3-aug-84 set error value
c               06-mar-86       removed return
c
c       ifail = 3
c     if(i.gt.30.and.wmx(i).gt.wlim) return
c
                wmx(i)=chop(wmx(i),wmin,wmax)
  800 continue
      ifail=0
      ero(6)=ers6                               !new
c
      if(kbug.eq.0) return
c
c               02-may-84       ibm to vax fortran change
c     call enkode('("after retrieval ..."/)',lbuf)
c     call enkode('("tskin =",f7.2/)',lbuf,tskin)
        write (6, 5030) tskin
5030    format (1x, 'after retrieval ...', /, 1x, 'tskin =', f7.2)
      do 820 k=1,nchan
                kc=k
                call vastau(tem,wmx,toto,vzen,tau,nsat,kc)
                call vasrte(tau,tem,tskin,rad,dbdt,tbc(k),dbdtbb,kc,is)
c                                               !changed lsta --> is
                dtb(k)=tbo(k)-tbc(k)
c
c               02-may-84       ibm to vax fortran change
c     call enkode('(" chan ",i2,", tbo =",f7.2,", tbc =",f7.2,", dtb =",
c    * f6.2/)',lbuf,kc,tbo(k),tbc(k),dtb(k))
        write(6, 5040) kc, tbo(k), tbc(k), dtb(k)
5040    format (1x, ' chan ', i2, ', tbo =', f7.2, ', tbc =', f7.2,
     1        ', dtb =', f6.2)
  820 continue
c
      return
      end
      subroutine solvex(gam,x,xit,xtx,xiv,nf,ne)
      implicit real*8 (a-h,o-z)
      real gam,x,xit
      dimension x(nf,ne),xit(nf,ne),xtx(nf,nf),xiv(nf,nf)
c   x-transpose-x
      do 140 j=1,nf
      do 140 i=1,nf
      sum=0.
      if(i.eq.j) sum=gam
      do 130 k=1,ne
  130 sum=sum+x(i,k)*x(j,k)
      xtx(i,j)=sum
  140 continue
c   inverse
      call symvrt(xtx,xiv,nf,nf)
c   inverse * x-transpose
      do 160 k=1,ne
      do 160 i=1,nf
      sum=0.
      do 150 j=1,nf
  150 sum=sum+xiv(i,j)*x(j,k)
      xit(i,k)=sum
  160 continue
      return
      end
      subroutine symvrt(a,s,kn,n)
c
c               02-may-84       sequence numbers removed
c
      real*8 a(kn,kn),s(kn,kn),sum
      s(1,1)=1.0/a(1,1)
      if(n-1)110,110,120
  110 return
  120 do 130 j=2,n
      s(1,j)=a(1,j)
  130 continue
      do 160 i=2,n
      im=i-1
      do 150 j=i,n
      sum=0.0
      do 140 l=1,im
      sum=sum+s(l,i)*s(l,j)*s(l,l)
  140 continue
      s(i,j)=a(i,j)-sum
  150 continue
      s(i,i)=1.0/s(i,i)
  160 continue
      do 170 i=2,n
      im=i-1
      do 170 j=1,im
  170 s(i,j)=0.0
      nm=n-1
      do 180 ii=1,nm
      im=n-ii
      i=im+1
      do 180 j=1,im
      sum=s(j,i)*s(j,j)
      do 180 k=i,n
      s(k,j)=s(k,j)-s(k,i)*sum
  180 continue
      do 210 j=2,n
      jm=j-1
      jp=j+1
      do 210 i=1,jm
      s(i,j)=s(j,i)
      if(jp-n)190,190,210
  190 do 200 k=jp,n
      s(i,j)=s(i,j)+s(k,i)*s(k,j)/s(k,k)
  200 continue
  210 continue
      do 250 i=1,n
      ip=i+1
      sum=s(i,i)
      if(ip-n)220,220,240
  220 do 230 k=ip,n
      sum=sum+s(k,i)*s(k,i)/s(k,k)
  230 continue
  240 s(i,i)=sum
  250 continue
      do 260 i=1,n
      do 260 j=i,n
      s(j,i)=s(i,j)
  260 continue
      return
      end
      subroutine vastau(temp,wmix,toto,zenang,tau,isat,kchan)
      logical co2,h2o,con,o3,gam,dry,newchn,newatm
      common/taucmp/co2,h2o,con,o3,gam,dry
      common/taucof/coef(500),flag(6)
      common/taumdl/tauco2(40),tauh2o(40),taucon(40),tauo3(40)
      dimension temp(40),wmix(40),tau(40),tsave(40),wsave(40)
      data tsave/40*0./,wsave/40*0./,ksave/-1/,isave/-1/,nlevs/40/
c
c               27-aug-84       call added because of linker/library
c                               problems
c
        call lblkdat
      newchn=isat.ne.isave.or.kchan.ne.ksave
      isave=isat
      ksave=kchan
      dt=0.
      dw=0.
      do 110 i=1,nlevs
      dt=dt+temp(i)-tsave(i)
      tsave(i)=temp(i)
      dw=dw+wmix(i)-wsave(i)
      wsave(i)=wmix(i)
      tauco2(i)=1.
      taucon(i)=1.
      tauh2o(i)=1.
  110 tauo3 (i)=1.
      theta=zenang*.0174533
      path=1./cos(theta)
      newatm=dt.ne.0..or.dw.ne.0.
      if(newchn) call pretav(isat,kchan)
      if(newatm) call preatv(temp,wmix)
      if(.not.dry) go to 120
      flag(2)=1.
      flag(3)=1.
      flag(5)=1.
  120 if(.not.co2) go to 130
      if(flag(4).ne.0.) go to 130
      call co2tav(path)
  130 if(.not.h2o) go to 140
      if(flag(5).ne.0.) go to 140
      call h2otav(path)
  140 if(.not.con) go to 150
      if(flag(1).ne.0..and.flag(2).ne.0..and.flag(3).ne.0.) go to 150
      call contav(wmix,path)
  150 if(.not.o3) go to 160
      if(flag(6).ne.0.) go to 160
      call o3tav(toto,path)
  160 continue
      do 170 i=1,nlevs
  170 tau(i)=tauco2(i)*taucon(i)*tauh2o(i)*tauo3(i)
      if(.not.gam) go to 180
      call gamtav(tau,kchan)
  180 continue
      return
      end
        subroutine lblkdat
c
c               27-aug-84       block data routine changed because
c                               of linker/library problems
c
c     block data
      logical co2,h2o,con,o3,gam,dry
      common/taucmp/co2,h2o,con,o3,gam,dry
      co2 = .true.
      h2o = .true.
      con = .true.
      o3 = .true.
      gam = .true.
      dry = .false.
c
c               27-aug-84       return added
c
        return
      end
      subroutine pretav(isat,kchan)
c
c                 1 7-apr-84    change in file name specification
c     real*8 file
c
      common/atmosf/pref(40),delp(40),pps(40),tref(40),
     * delt(40),delts(40),deltss(40),tempf(40),tts(40),wvap(40),nl,nm,nw
      dimension px(40),tx(40),cbuf(506,12)
      common/taucof/buf(506)
      data px/.1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,50.,60.
     *,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,430.,475.,
     * 500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./
      data tx/235.51,250.04,266.99,270.65,265.04,257.88,249.45,243.64,
     * 239.23,232.73,227.70,225.02,223.13,221.68,220.50,217.23,216.83,
     * 7*216.65,220.79,228.58,235.47,241.44,244.83,249.51,251.92,258.36,
     * 262.49,266.37,268.62,274.18,278.73,282.96,284.68,287.43/
      data iu/13/,init/0/,mrec/13/
c
c            1 7-apr-84 change in file name specification
c     data file/'vastaucf'/,len/506/
c
        character*200 file
        integer len_file
c        data file /'../static/goeslib/for044.dat'/
c     byte file(28)
c     data file/'',
c     1   'v','a','s','t','a','u','c','f','.','d','a','t',"0/
      data len/506/
      call get_directory('static',file,len_file)
      file = file(1:len_file)//'goeslib/for044.dat'

      if(init.ne.0) go to 150
      nl=40
      nw=23
      nm=nl-1
      nw=nl-nw+1
      do 110 i=1,nl
      pref(i)=px(i)
  110 tref(i)=tx(i)
      delp(1)=pref(1)
      do 120 i=1,nm
  120 delp(i+1)=pref(i+1)-pref(i)
      pps(1)=0.
      do 140 i=1,nm
  140 pps(i+1)=alog((pref(i+1)+pref(i))/2000.)
  150 if(isat.eq.init) go to 200
      call dopen(file,iu,len)
      irec=(isat-1)*mrec
      nrec=mrec-1
      do 160 j=1,nrec
      irec=irec+1
      call dread(iu,irec,cbuf(1,j))
  160 continue
      call dclose(iu)
      init=isat
  200 do 210 i=1,len
  210 buf(i)=cbuf(i,kchan)
      return
      end
      subroutine preatv(temp,wmix)
      common/atmosf/pref(40),delp(40),pps(40),tref(40),
     * delt(40),delts(40),deltss(40),tempf(40),tts(40),wvap(40),nl,nm,nw
      dimension temp(40),wmix(40)
      data cgrav/.509858e-3/
      nlw=nw-1
      sum1=0.
      sum2=0.
      do 110 i=1,nl
      delt(i)=temp(i)-tref(i)
      sum1=sum1+delp(i)*delt(i)
      sum2=sum2+delp(i)*delt(i)*pref(i)
      delts(i)=sum1/pref(i)
      deltss(i)=2.*sum2/(pref(i)*pref(i))
  110 tempf(i)=1./temp(i)
      tts(1)=0.
      do 120 i=1,nm
  120 tts(i+1)=alog((temp(i)+temp(i+1))/546.)
      call ulmr(wmix)
      do 130 i=1,nlw
  130 wvap(i)=0.
      wvap(nw)=cgrav*delp(nw)*wmix(nw)
      do 140 i=nw,nm
  140 wvap(i+1)=wvap(i)+cgrav*delp(i+1)*(wmix(i)+wmix(i+1))
      return
      end
      subroutine ulmr(w)
      dimension w(40),pw(8),wp(8)
      data pw/70.,85.,100.,115.,135.,150.,200.,250./
      data init/1/
      if(init.eq.0) go to 110
      do 100 i=1,8
      wp(i)=(pw(i)/300.)**3
  100 continue
      init=0
  110 wb=w(26)
      do 120 i=1,8
      j=i+17
      w(j)=wb*wp(i)
  120 continue
      return
      end
      subroutine co2tav(path)
      common/taucof/co2cof(5,40),co2ang(3,40),concof(6),h2ocof(14)
     *,o3cof(2,40),o3ang(2,40),flag(6)
      common/atmosf/pref(40),delp(40),pps(40),tref(40),
     * delt(40),delts(40),deltss(40),tempf(40),tts(40),wvap(40),nl,nm,nw
      common/taumdl/tauco2(40),tauh2o(40),taucon(40),tauo3(40)
      sec=path-1.
      secs=sec*sec
      tau=1.
      do 110 j=1,nl
      tau=tau*(co2cof(1,j)*delt(j)+co2cof(2,j)*delt(j)*delt(j)
     * +co2cof(3,j)*delts(j)+co2cof(4,j)*deltss(j)+co2cof(5,j))
c
c            1 7-apr-84 correction of transmission error
c     tauzu
c
        tauz = tau
      if(sec.eq.0.) go to 110
      tauz=tauz+(co2ang(1,j)*sec+co2ang(2,j)*sec*deltss(j)
     *      +co2ang(3,j)*secs)
  110 tauco2(j)=tauz
      return
      end
      subroutine h2otav(path)
      common/taucof/co2cof(5,40),co2ang(3,40),concof(6),h2ocof(14)
     *,o3cof(2,40),o3ang(2,40),flag(6)
      common/atmosf/pref(40),delp(40),pps(40),tref(40),
     * delt(40),delts(40),deltss(40),tempf(40),tts(40),wvap(40),nl,nm,nw
      common/taumdl/tauco2(40),tauh2o(40),taucon(40),tauo3(40)
      dimension pvap(40)
      data pvap/40*0./
      data thresh/1.e-4/
      do 110 k=nw,nl
  110 pvap(k)=wvap(k)*path
      iflag=0
      do 240 k=nw,nl
      if(pvap(k).lt.thresh) go to 240
      a3=pps(k)
      a4=tts(k)
      a9=a3*a4
      a12=a4*a4
      cof0=h2ocof(1)+a4*h2ocof(4)+a3*h2ocof(3)+a12*h2ocof(12)
     * +a9*h2ocof(9)
      cof1=h2ocof(2)+a4*h2ocof(6)+a3*h2ocof(5)+a12*h2ocof(11)
     * +a9*h2ocof(13)
      cof2=h2ocof(7)+a4*h2ocof(8)+a3*h2ocof(14)
      cof3=h2ocof(10)
      iflag=iflag+1
      if(iflag.ne.1)go to 120
      ustart=pvap(k)
      go to 170
  120 continue
      a22=a2*a2
      a23=a2*a22
      wsum2=cof0+cof1*a2+cof2*a22+cof3*a23
      wsum1=wsum2-wsum
      deriv=cof1+2.0*cof2*a2+3.0*cof3*a22
      if(deriv.gt.1.) go to 130
      ustart=pvap(k)
      go to 170
  130 continue
      if(abs(deriv).gt..0001)go to 140
      go to 150
  140 a2=a2-wsum1/deriv
  150 arg=10.*a2-a4
      d=0.
      if(arg.lt.-100.) go to 160
      d=exp(arg)
  160 ustart=d+pvap(k)-pvap(k-1)
  170 if(ustart.gt.thresh) go to 180
      ustart=pvap(k)
  180 a2=.1*(alog(ustart)+a4)
      a22=a2*a2
      a23=a2*a22
      wsum=cof0+cof1*a2+cof2*a22+cof3*a23
      if(wsum.gt.5.) go to 190
      if(wsum.lt.-20.) go to 200
      go to 210
  190 continue
      tau=0.0
      go to 220
  200 continue
      tau=1.0
      go to 220
  210 continue
      tau=exp(-exp(wsum))
  220 continue
      if(tau.gt.1.0)go to 240
      if(tau.gt.0.0)go to 230
      tau=0.0
  230 tauh2o(k)=tau
  240 continue
      return
      end
      subroutine contav(wmix,path)
      common/taucof/co2cof(5,40),co2ang(3,40),concof(6),h2ocof(14)
     *,o3cof(2,40),o3ang(2,40),flag(6)
      common/atmosf/pref(40),delp(40),pps(40),tref(40),
     * delt(40),delts(40),deltss(40),tempf(40),tts(40),wvap(40),nl,nm,nw
      common/taumdl/tauco2(40),tauh2o(40),taucon(40),tauo3(40)
      dimension trn(40),wmix(40)
      do 110 i=1,nl
  110 trn(i)=0.
      if(flag(1).ne.0.) go to 130
      tau=0.
      f1=0.
      do 120 i=2,nl
      f2=pref(i)*tempf(i)
      tau=tau+(f1+f2)*delp(i)/2.
      f1=f2
  120 trn(i)=tau*concof(1)
  130 if(flag(2).ne.0.) go to 150
      tau=0.
      f1=0.
      do 140 i=nw,nl
      f2=pref(i)*wmix(i)
      tau=tau+(f1+f2)*delp(i)/2.
      f1=f2
  140 trn(i)=trn(i)+tau*concof(3)
  150 if(flag(3).ne.0.) go to 170
      tau=0.
      f1=0.
      do 160 i=nw,nl
      x=concof(6)*tempf(i)
      f2=exp(x)*pref(i)*wmix(i)**2
      tau=tau+(f1+f2)*delp(i)/2.
      f1=f2
  160 trn(i)=trn(i)+tau*concof(5)
  170 do 180 i=2,nl
      if(trn(i).eq.0.) go to 180
      taucon(i)=exp(-path*trn(i))
  180 continue
      return
      end
      subroutine o3tav(toto,path)
      common/taucof/co2cof(5,40),co2ang(3,40),concof(6),h2ocof(14)
     *,o3cof(2,40),o3ang(2,40),flag(6)
      common/taumdl/tauco2(40),tauh2o(40),taucon(40),tauo3(40)
      data nl/40/
      sec=path-1.
      frac2=(toto-257.)/223.
      frac1=1.-frac2
      do 110 j=1,nl
      arg1=o3cof(1,j)
      arg2=o3cof(2,j)
      if(sec.eq.0.) go to 110
      exp1=1.+sec*o3ang(1,j)
      exp2=1.+sec*o3ang(2,j)
      arg1=arg1**exp1
      arg2=arg2**exp2
  110 tauo3(j)=frac1*arg1+frac2*arg2
      return
      end
      subroutine gamtav(tau,kchan)
      dimension tau(*)
      common/gde/gv(12),dv(12),ev(12)
      data taumin/.0001/,taumax/.9999/,nl/40/
      g=gv(kchan)
      if(g.eq.0.) return
      if(g.eq.1.) return
      do 100 i=1,nl
      if(tau(i).le.taumin) go to 100
      if(tau(i).ge.taumax) go to 100
      tau(i)=tau(i)**g
  100 continue
      return
      end
      subroutine plnkiv(isat)
c
c            1 7-apr-84 change in file name specification
c     real*8 file
c
      common/plankv/fnu(12),fk1(12),fk2(12),tc(2,12)
      common/gde/gv(12),dv(12),ev(12)
      common/tsurfc/tsc(4,2)
      common/tsurfe/tse(4,2)
      common/use/iuch(12,2)
      common/chkcof/chkc(5)
      dimension buf(506)
c
c            1 7-apr-84 change in file name specification
c     data file/'vastaucf'/,len/506/
c
        character*200 file
        integer len
c        data file /'../static/goeslib/for044.dat'/
c      byte file(28)
c      data file/'',
c     1   'v','a','s','t','a','u','c','f','.','d','a','t',"0/
      data len/506/
      data iu/13/,nc/12/
c
c            1 3-jun-j4 debug
c
c               type 1090, isat

      write(6,*)' PLNKIV 1'
      call get_directory('static',file,len_file)
      file = file(1:len_file)//'goeslib/for044.dat'
      call s_len(file,len_file)
1090            format (/, 1x, ' satellite ', i3)
      call dopen(file(1:len_file),iu,len)
      irec=isat*(nc+1)
      call dread(iu,irec,buf)
      call dclose(iu)
      mc=nc*2
      do 110 j=1,nc
      fnu(j)=buf(j)
      fk1(j)=buf(j+nc)
  110 fk2(j)=buf(j+mc)
c
c            1 3-jun-84 debug
c
c               type 1000, (fnu(iii), fk1(iii), fk2(iii), iii = 1, nc)
1000            format (/, 1x, 'fnu         fk1          fk2', /,
     1                (3e12.4))
      n=nc*3
      write(6,*)' PLNKIV 2'
      do 130 j=1,nc
      do 120 i=1,2
      n=n+1
  120 tc(i,j)=buf(n)
  130 continue
      write(6,*)' PLNKIV 3'
c
c            1 3-jun-84 debug
c
c               type 1010, ((tc(iii, jjj), iii = 1, 2), jjj = 1, nc)
1010            format (/, 1x, 'tc1        tc2', /, (2e12.4))
      do 140 i=1,nc
      n=n+1
  140 gv(i)=buf(n)
      write(6,*)' PLNKIV 3a'
c
c            1 3-jun-84 debug
c
c               type 1020, (gv(iii), iii = 1, nc)
1020            format (/, 1x, 'gb', /, (e12.4))
      do 150 i=1,nc
      n=n+1
  150 dv(i)=buf(n)
      write(6,*)' PLNKIV 3b'
c
c            1 3-jun-84 debug
c
c               type 1030, (dv(iii), iii = 1, nc)
1030            format (/, 1x, 'dv', /, (e12.4))
      do 160 i=1,nc
      n=n+1
  160 ev(i)=buf(n)
      write(6,*)' PLNKIV 3c'
c
c            1 3-jun-84 debug
c
c               type 1040, (ev(iii), iii = 1, nc)
1040            format (/, 1x, 'ev', /, (e12.4))
      n=100
      do 180 j=1,2
      do 170 i=1,4
      n=n+1
  170 tsc(i,j)=buf(n)
  180 continue
      write(6,*)' PLNKIV 4'
c
c            1 3-jun-84 debug
c
c               type 1050, ((tsc(iii, jjj), iii = 1, 2), jjj = 1, 4)
1050            format (/, 1x, 'tsc1          tsc2', /,
     1                (2e12.4))
      n=110
      do 200 j=1,2
      do 190 i=1,nc
      n=n+1
  190 iuch(i,j)=buf(n)
  200 continue
c
c            1 3-jun-84 debug
c
c               type 1060, ((iuch(iii, jjj), jjj = 1, 2), iii = 1, nc)
1060            format (/, 1x, 'iuch1         iuch2', /,
     1                (2i12))
      n=140
      do 210 i=1,5
      n=n+1
  210 chkc(i)=buf(n)
      n=150
c
c            1 3-jun-84 deubg
c
c               type 1070, (chkc(iii), iii = 1, 5)
1070            format (/, 1x, 'chkc', /, (5e12.4))
      do 230 j=1,2
      do 220 i=1,4
      n=n+1
  220 tse(i,j)=buf(n)
  230 continue
      write(6,*)' PLNKIV 5'
c
c            1 3-jun-84 debug
c
c               type 1080, ((tse(iii, jjj), jjj = 1, 2), iii = 1, 4)
1080            format (/, 1x, 'tse1         tse2', /,
     1                (2e12.4))
      return
      end
      function vplanc(t,k)
      common/plankv/fnu(12),fk1(12),fk2(12),tc(2,12)

      
      tt=tc(1,k)+tc(2,k)*t
      expn=exp(fk2(k)/tt) - 1.
      vplanc=fk1(k)/expn
      return
      end
      function vbrite(r,k)
      common/plankv/fnu(12),fk1(12),fk2(12),tc(2,12)
      expn=fk1(k)/r+1.
      tt=fk2(k)/alog(expn)
      vbrite=(tt-tc(1,k))/tc(2,k)
      return
      end
      function vdbdtb(tbb,k)
      common/plankv/fnu(12),fk1(12),fk2(12),tc(2,12)
      t=tbb
      tt=tc(1,k)+tc(2,k)*t
      expn=exp(fk2(k)/tt) - 1.
      b=fk1(k)/expn
      bb=b*b
      tt=t*t
      expn=exp(fk2(k)/t)
      q=fk2(k)/fk1(k)
      vdbdtb=q*expn*bb/tt
      return
      end
      subroutine vasrte(tau,temp,tsfc,rad,dbdt,tbb,dbdtbb,kchan,lgnd)
      common/plankv/fnu(12),fk1(12),fk2(12),tc(2,12)
      common/gde/gv(12),dv(12),ev(12)
      dimension tau(*),temp(*),dbdt(*)
      data nl/40/
      f1=fk1(kchan)
      f2=fk2(kchan)
      q=f2/f1
      t1=temp(1)
      b1=vplanc(t1,kchan)
      bb=b1*b1
      tt=t1*t1
      ex=q*exp(f2/t1)
      dbdt(1)=ex*bb/tt
      tau1=tau(1)
      rad=0.
      do 110 i=2,nl
      t2=temp(i)
      b2=vplanc(t2,kchan)
      bb=b2*b2
      tt=t2*t2
      ex=q*exp(f2/t2)
      dbdt(i)=ex*bb/tt
      tau2=tau(i)
      if(i.gt.lgnd) go to 100
      rad=rad+.5*(b1+b2)*(tau1-tau2)
c
c               02-may-84       debug
c
c               type 1000, i, b1, b2, tau1, tau2, rad
1000            format (i3, 5(1x, e12.5))
  100 b1=b2
      tau1=tau2
  110 continue
      bs=vplanc(tsfc,kchan)
      rad=rad+bs*tau(lgnd)
      rad=rad+dv(kchan)
      rad=amax1(rad,.001)
      tbb=vbrite(rad,kchan)
      tt=tbb*tbb
      bb=rad*rad
      ex=q*exp(f2/tbb)
      dbdtbb=ex*bb/tt
      return
      end
        function chop(c, a, b)
c
c               04-may-84       function added
c               06-jun-84       function modified
c
        x = amax1(c, a)
        x = amin1(x, b)
        chop = x
        return
        end
        subroutine dclose (iu)
c
c            1 7-apr-84 subroutine added
c
        common/dio/lenn
        close (unit=iu)
        return
        end
        subroutine dopen (file, iu, len)
c
c            1 7-apr-84 subroutine added
c
        common/dio/lenn
c       byte file(1)
        character*200 file
c
        lenn = len
        
        open (unit=iu, file=file,
     1      form='formatted')
c
        return
        end
        subroutine dread (iu, irec, buf)
c
c            1 7-apr-84 subroutine added
c
        common/dio/lenn
        dimension buf(lenn), allbuf(506,39), global (19734)
        equivalence (allbuf(1,1),global(1))
        integer test
        data test /0/
c
c       read (unit=iu, rec=irec) (buf(i), i = 1, lenn)
        if (test.eq.1) then
        continue
        else
        test = 1
        read (iu, *) (global(i), i=1,19734)
        endif

c       go here directly if the file has been read before
        do i = 1,lenn
        buf(i) = allbuf(i,irec)
        enddo

c
        return
        end
        subroutine wmix (p, t, dd, w, nl)
c
c               07-jun-84       source code from wang
c                               (not on u. of w. source tape)
c
c               dd      - dewpoint depression
c     modified dimension structure for sun machine warning error 1/16/02 DB
      
        integer nl
        dimension p(nl), t(nl), dd(nl), w(nl)
        do 130 i = 1, nl
           td = t(i) - dd(i)
           if (td .gt. 253.0) go to 110
           es = vpice(td)
           go to 120
 110       es = satvap(td)
 120       w(i) = 622.0 * es / p(i)
 130    continue
        return
        end

      function wsat(p, t)
c     
c     1 4-may-84 added from csu vas code
c     (not on u. of w. source tape)
c     modified by adding dd(1) to make compliant with revised
c     function parameter list 1/16/02 db
      real dd(1)
      dd(1) = 0.0
c     
      call wmix(p,t,dd,w,1)
      wsat=w
c     
c     1 4-may-84 return added
c     
      return
      end
      
      function vpice(temp)
c
c               07-jun-84       source code from wong
c                               (not on u. of w. source tape)
c
        dimension c(6)
c
        real*8 c, t
c
        data c/ 0.7859063157d+00,  0.3579242320d-01,
     1       -0.1292820828d-03,  0.5937519208d-06,
     1        0.4482949133d-09,  0.2176664827d-10/
c
        t = temp - 273.16
        vplog = c(1) + t * (c(2) + t * (c(3) + t * (c(4) + t * (c(5) +
     1               t *  c(6)))))
        vpice = 10.0 ** vplog
        return
        end
        function satvap(temp)
c
c               07-jun-84       source code from wang
c                               (not on u. of w. source tape)
c
        dimension c(10)
c
        real*8 c, s, t
c
        data c/ 0.9999968760d-00, -0.9082695004d-02,
     1        0.7873616869d-04, -0.6111795727d-06,
     1        0.4388418740d-08, -0.2988388486d-10,
     1        0.2187442495d-12, -0.1789232111d-14,
     1        0.1111201803d-16, -0.3099457145d-19/
        t = temp - 273.16
        s = 6.1078000000d+00 / (c(1) + t * (c(2) + t * (c(3) +
     1                               t * (c(4) + t * (c(5) +
     1                               t * (c(6) + t * (c(7) +
     1                               t * (c(8) + t * (c(9) +
     1                               t * c(10)))))))))) ** 8
        satvap = s
        return
        end
      subroutine precw(p,w,u,np)
c
c               02-may-84       sequence numbers removed
c
      dimension p(*),w(*),u(*)
      data f/1961.33/
      w1=w(1)
      p1=p(1)
      s=0.
      u(1)=s
c
c               07-jun-84       correction of tape read error
c     do 100 i=2,n
        do 100 i = 2, np
      w2=w(i)
      p2=p(i)
      dp=abs(p2-p1)
      s=s+(w1+w2)*dp/f
      u(i)=s
      w1=w2
  100 p1=p2
      return
      end
      function vskint(tbb,jcof,juse,jact)
c * obtain sfc-skin-temp from vas brightness temperatures
c        jcof=0 for empirical coefficients
c        jcof=1 for climatological(theoretical) coeffs.
c
c        juse=0 to use 'day'(2-chan) or 'night'(3-chan) eqn.
c           according to solar zenith angle
c *******  currently forced out; will use day only if chan 12 missing.
c        juse=1 to use 'day' eqn.
c        juse=2 to use 'night' eqn.
c
c        jact=0 if nothing done
c        jact=1 if 'day' eqn. actually used
c        jact=2 if 'night' eqn. actually used
      common/nav/vlat,vlon,vzen,szen,il,ie,iras,ipic,itime,jtime,jday
      common/tsurfc/tscc(4,2)
      common/tsurfe/tsce(4,2)
      dimension tsc(4,2),tx(3),ic(3),tbb(*)
      data ic/7,8,12/,nx/3/,vmisg/999999./
      do 90 i=1,nx
      j=ic(i)
   90 tx(i)=tbb(j)
      jact=0
      ts=tx(2)
      if(tx(1).eq.vmisg) go to 180
      if(tx(2).eq.vmisg) go to 180
      if(jcof.ne.0) go to 110
      do 100 j=1,2
      do 100 i=1,4
  100 tsc(i,j)=tsce(i,j)
      go to 130
  110 do 120 j=1,2
      do 120 i=1,4
  120 tsc(i,j)=tscc(i,j)
  130 j=juse
      if(j.ne.0) go to 150
      j=2                       !changed 1 --> 2
c     if(szen.ge.90.) j=2       !deleted line
  150 if(j.eq.1) go to 160
      if(tx(3).eq.vmisg) j=1
      if(tx(3).lt.tx(2)-4.) j=1 !changed tx(2) --> tx(2)-4.
  160 jact=j
      ts=tsc(4,j)
      do 170 i=1,nx
  170 ts=ts+tx(i)*tsc(i,j)
      if(ts.eq.0.) ts=vmisg
  180 vskint=ts
      return
      end
