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

      subroutine taugim(t,w,o,theta,ngoes,kan,tau)
c * GOES/I-M TAU via McMillin-Fleming-Woolf-Eyre regression model
c   Input temperatures ('T'), water-vapor mixing ratios ('W'),
c        and ozone mixing ratios ('O') must be at the 'STANDARD'
c        40 levels used for radiative-transfer calculations.
c   If ozone profile is not available, pass a zero-filled array.  This
c        will cause the reference profile to be used.
c   Units -- T: degrees Kelvin; W: grams/kilogram; O: ppm by volume.
c   Logical unit number 15 is used for the coefficient file.
c      THETA = local (satellite) zenith angle in degrees
c      NGOES = GOES satellite number, e.g. 8 (GOES/I)
c        KAN = channel number (1,...,25)
c               1 thru 18 are SOUNDER
c              22 thru 25 are IMAGER
c   To get Planck-function, band correction, and TSKIN coefficients,
c    plus gammas, deltas, epsilons, and 'use' flags,
c        user *MUST* call PFCGIM!
c   Total transmittance is returned via the formal parameter TAU.
c   Separate DRY/WET/OZO transmittances are returned via COMMON.
c
      parameter (nc=3,nk=25,nl=40,nlm1=nl-1,nt=2,nv=20,nx=9,nxp1=nx+1)
      parameter (iuc=15,lenc=nxp1*nl*nc,lencb=lenc*4)
      common/refatm/pref(nl),tref(nl),wref(nl),oref(nl)
      common/taudwo/taud(nl),tauw(nl),tauo(nl)
      dimension t(nl),w(nl),o(nl),tau(nl),tl(nl),wl(nl),ol(nl)
      dimension delo(nl),delt(nl),delw(nl)
      dimension sumo(nl),sums(nl),sumt(nl),sumw(nl)
      dimension dp(nl),prsq(nl),xx(nx,nl,nc)
      dimension coef(nxp1,nl,nc,nk),tauc(nl,nc),sqp(nc)
      dimension wdep(nl),sqw(nl),odep(nl),sqo(nl),sqd(nl),sqdep(nl,nc)
      equivalence (sqd(1),sqdep(1,1)),(taud(1),tauc(1,1))
      equivalence (sqw(1),sqdep(1,2)),(tauw(1),tauc(1,2))
      equivalence (sqo(1),sqdep(1,3)),(tauo(1),tauc(1,3))
      save
      real fast_coef(30000)
      equivalence (fast_coef(1), coef(1,1,1,1) )
      character*8 cfile
      data cfile/'GOESRTCF'/
      character*200 fname
      integer len
c
      logical oldatm,oldang
      data tl/nl*0./,wl/nl*0./,ol/nl*0./,init/0/
      data sqd/nl*1./,trap/-999.99/
      secant(z) = 1./cos(0.01745329*z)
c
      if(kan.gt.18.and.kan.lt.22) go to 200
      if(init.eq.ngoes) go to 100



      call get_directory('static',fname,len)
      if(ngoes.eq.8) then
         open(23,file=fname(1:len)//'fastcoef8.dat',form='formatted')
      elseif(ngoes.eq.9) then
         open(23,file=fname(1:len)//'fastcoef9.dat',form='formatted')
      endif
      read(23,*) fast_coef
      close(23)



c      open(iuc,file=cfile,recl=lencb,access='direct',status='old')
c      irec=(ngoes-8)*nk
c      do l=1,nk
c         irec=irec+1
c         read(iuc,rec=irec) (((coef(i,j,k,l),i=1,nxp1),j=1,nl),k=1,nc)
c      enddo
c      close(iuc)
      prsq(1)=pref(1)**2
      dp(1)=pref(1)
      do j=1,nlm1
         jp1=j+1
         prsq(jp1)=pref(jp1)**2
         dp(jp1)=pref(jp1)-pref(j)
      enddo
      thetl=trap
      init=ngoes
c
  100 dt=0.
      dw=0.
      do=0.
      if(o(1).eq.0.) then
         do j=1,nl
            o(j)=oref(j)
         enddo
      endif
      do j=1,nl
crma     dt=dt+abs(t(j)-tl(j))
         dt=dt+t(j)-tl(j)
         tl(j)=t(j)
crma     dw=dw+abs(w(j)-wl(j))
         dw=dw+w(j)-wl(j)
         wl(j)=w(j)
crma     do=do+abs(o(j)-ol(j))
         do=do+o(j)-ol(j)
         ol(j)=o(j)
      enddo
      dtwo=dt+dw+do
      oldatm=dtwo.eq.0.
      if(oldatm) go to 110
      st=0.
      ss=0.
      sw=0.
      so=0.
      do j=1,nl
         delt(j)=t(j)-tref(j)
         delw(j)=w(j)-wref(j)
         delo(j)=o(j)-oref(j)
         odep(j)=dp(j)*o(j)
         wdep(j)=dp(j)*w(j)
         if(j.gt.1) then
            jm1=j-1
            delt(j)=0.5*(delt(j)+t(jm1)-tref(jm1))
            delw(j)=0.5*(delw(j)+w(jm1)-wref(jm1))
            delo(j)=0.5*(delo(j)+o(jm1)-oref(jm1))
            odep(j)=dp(j)*0.5*(o(j)+o(jm1))
            if(j.gt.nv) then
               wdep(j)=dp(j)*0.5*(w(j)+w(jm1))
            else
               wdep(j)=wdep(jm1)+dp(j)*0.5*(w(j)+w(jm1))
            endif
         endif
         sqo(j)=sqrt(odep(j))
         sqw(j)=sqrt(wdep(j))
         dtdp=delt(j)*dp(j)
         st=st+dtdp
         sumt(j)=st/pref(j)
         ss=ss+dtdp*pref(j)
         sums(j)=2.*ss/prsq(j)
         dwdp=delw(j)*dp(j)
         sw=sw+dwdp*pref(j)
         sumw(j)=2.*sw/prsq(j)
         dodp=delo(j)*dp(j)
         so=so+dodp*pref(j)
         sumo(j)=2.*so/prsq(j)
      enddo
  110 oldang=theta.eq.thetl
      if(.not.oldang) then
         if(theta.eq.0.) then
            path=1.
            sqrtp=1.
         else
            path=secant(theta)
            sqrtp=sqrt(path)
         endif
         pathm1=path-1.
         sqp(1)=1.
         sqp(2)=sqrtp
         sqp(3)=sqrtp
         thetl=theta
      endif
      if(oldatm.and.oldang) go to 120
      do j=1,nl
         dt=delt(j)
         ss=sums(j)
         st=sumt(j)
c **** DRY
         xx(1,j,1)=dt*path
         xx(2,j,1)=dt*dt*path
         xx(3,j,1)=st*path
         xx(4,j,1)=ss*path
         xx(5,j,1)=pathm1
         xx(6,j,1)=pathm1*pathm1
         xx(7,j,1)=st*pathm1
         xx(8,j,1)=ss*pathm1
         xx(9,j,1)=dt*pathm1
c **** WET
         sqpw=sqrtp*sqw(j)
         dw=delw(j)
         sw=sumw(j)
         xx(1,j,2)=dt
         xx(2,j,2)=ss
         xx(3,j,2)=dw
         xx(4,j,2)=sw
         xx(5,j,2)=dt*sqpw
         xx(6,j,2)=dt*dt*sqpw
         xx(7,j,2)=dw*sqpw
         xx(8,j,2)=dw*dw*sqpw
         xx(9,j,2)=dw*dt*sqpw
c **** OZO
         sqpo=sqrtp*sqo(j)
         do=delo(j)
         so=sumo(j)
         xx(1,j,3)=dt
         xx(2,j,3)=ss
         xx(3,j,3)=do
         xx(4,j,3)=so
         xx(5,j,3)=dt*sqpo
         xx(6,j,3)=dt*dt*sqpo
         xx(7,j,3)=do*sqpo
         xx(8,j,3)=do*do*sqpo
         xx(9,j,3)=do*dt*sqpo
      enddo
  120 do k=1,nc
         do j=1,nl
            tauc(j,k)=1.0
         enddo
      enddo
c
      l=kan
      do k=1,nc
         taul=1.
         do 130 j=1,nl
         if(taul.eq.0.) go to 130
         yy=coef(nxp1,j,k,l)
         if(yy.eq.trap) then
            taul=0.
            go to 130
         endif
         do i=1,nx
            yy=yy+coef(i,j,k,l)*xx(i,j,k)
         enddo
         yy=amin1(yy,0.)
         tauy=taul*exp(yy*sqdep(j,k)*sqp(k))
         taul=tauy
  130    tauc(j,k)=taul
      enddo
      do j=1,nl
         tau(j)=taud(j)*tauw(j)*tauo(j)
      enddo
  200 return
      end
      block data refpro
c $ Transmittance-Model Reference = U.S. Standard Atmosphere, 1976
      parameter (nl=40)
      common/refatm/pref(nl),tref(nl),wref(nl),oref(nl)
      data pref/ .1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,
     + 50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
     + 430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./
      data tref/
     +  231.70, 245.22, 263.35, 270.63, 264.07, 257.93, 249.51, 243.65,
     +  239.24, 232.64, 228.07, 225.00, 223.13, 221.72, 220.54, 217.28,
     +  216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.72,
     +  220.85, 228.58, 235.38, 241.45, 244.81, 249.48, 251.95, 258.32,
     +  262.48, 266.40, 268.61, 274.21, 278.74, 282.97, 284.71, 287.50/
      data wref/
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.004,  0.005,  0.014,
     +   0.036,  0.089,  0.212,  0.331,  0.427,  0.588,  0.699,  1.059,
     +   1.368,  1.752,  1.969,  2.741,  3.366,  3.976,  4.255,  4.701/
c    +   0.004,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,
c    +   0.005,  0.005,  0.005,  0.005,  0.005,  0.004,  0.004,  0.004,
c    +   0.004,  0.004,  0.004,  0.004,  0.005,  0.006,  0.008,  0.022,
c    +   0.057,  0.143,  0.340,  0.533,  0.687,  0.946,  1.125,  1.704,
c    +   2.202,  2.819,  3.168,  4.411,  5.416,  6.397,  6.846,  7.564/
      data oref/
     + 0.65318,1.04797,2.13548,3.82386,5.26768,6.11313,7.35964,7.75004,
     + 7.82119,7.56126,6.92006,6.10266,5.55513,5.15298,4.59906,2.86792,
     + 2.29259,1.80627,1.28988,0.93973,0.72277,0.54848,0.46009,0.29116,
     + 0.16277,0.09861,0.06369,0.05193,0.04718,0.04097,0.03966,0.03614,
     + 0.03384,0.03342,0.03319,0.03249,0.03070,0.02878,0.02805,0.02689/
      end
      subroutine pfcgim(ngoes)
c * Input GOES/I-M Planck-function, band-cor'n & TSKIN coeff's
c        plus gammas, deltas, epsilons, and 'use' flags
c        NGOES = GOES satellite number, e.g. 8 (GOES/I)
      parameter (iuc=15,lenc=1200,nk=25,nt=2,lenp=nk*(nt+3),lent=30)
      parameter (leng=nk*3,lenu=nk*nt,lenr=lenc*4)
      parameter (mbg=200,mbu=300,nrc=20)
      common/gimgde/gbuf(leng)
      common/plncgo/pbuf(lenp)
      common/tskcof/tbuf(lent)
      common/use/ibuf(lenu)
      dimension cbuf(lenc)
      save
      character*8 cfile
      data cfile/'GOESRTCF'/
      character*200 fname
      integer len
c
c
c
c
c
c....begin portable code
c
c
      call get_directory('static',fname,len)
      if(ngoes.eq.8) then
         open(24,file=fname(1:len)//'pfcgim8.dat',form='formatted')
      elseif(ngoes.eq.9) then
         open(24,file=fname(1:len)//'pfcgim9.dat',form='formatted')
      endif
      read(24,*) pbuf,tbuf,gbuf,ibuf
      close(24)

c
c
c....dont call original code (follows if test)... only runs on IBM
c
c
      if (1.eq.1) return
c
c

      open(iuc,file=cfile,recl=lenr,access='direct',status='old')
      irc=(ngoes-8)*nk+nrc
      read(iuc,rec=irc) cbuf
      close(iuc)
      do l=1,lenp
         pbuf(l)=cbuf(l)
      enddo
      m=lenp
      do l=1,lent
         m=m+1
         tbuf(l)=cbuf(m)
      enddo
      m=mbg
      do l=1,leng
         m=m+1
         gbuf(l)=cbuf(m)
      enddo
      m=mbu
      do l=1,lenu
         m=m+1
         ibuf(l)=cbuf(m)
      enddo
      return
      end
