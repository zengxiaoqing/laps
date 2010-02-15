
module src_versuch



 use src_global
 use src_hilfe
 use src_hyd
 use src_scatter
 use src_tmatrix
 use MODI_GAMMA

!******************************************************************
!Monika Pfeifer 28.09.2004
!INTENT(IN) :: rain, snow, graup
!INTENT(IN) :: temperatur, druck, qd
!INTENT(OUT) :: Zhhx, Ahhx, ZDR, LDR
!******************************************************************

contains
subroutine versuch(rain,snow,graup,tmp,zhhx,ahhx,zhvx, zvhx, zvvx, druck, qd, delahvx, kkdp, rhohvx)

Implicit none


!Include 'variablen.incf'

!common/glramssd/glszd
!common/glhydbdy/glbdy,glcad
!common/cantpara/cdpnum,cdprb,cdpth,&
!          cdpab,cdpdd,cdpab2,cdpdd2,cdpabd

REAL :: cloud, rain, snow, graup, tmp,  zhhx, ahhx,zzzdr, llldr, druck, qd, &
     thetag, thetagd, thetas, thetasd, delahvx, kkdp, rhohvx, &!
     zvhx, zhvx, zvvx
REAL :: qh, xx, erfct, titu, crtnm, wvln, qitu, wt(40), qg  !

REAL :: arg1, arg2

Complex  :: epswtr, epsice, epsplt
character :: sctmodel*3
integer, parameter :: msz=40

Integer  :: ii, idr, nrank,&
            n0rank, n0pnt, ist, nsz, isz, iisz !

Real :: sz(msz),sdprb(msz), thrdr(2)

Real :: freq0, lambd(9), tv, qs, scmixgrice, t, q, epsreplt, epsimplt !

Integer :: ssznum  !
Real*4 :: temp, rhox, testgm !
CHARACTER*9 :: hyddo !



! equivalence (glbdy(2),temp),(glbdy(3),scmix),&
!       (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp),&
!       (glbdy(7),deq),(glbdy(8),dmx),&
!       (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)

integer :: &
nn

real :: &
n0s, &
n0s1, &
n0s2, &
alf, &
bet, &
hlp, &
zams, &
m2s, &
m3s, &
tm, &
lambda

real, dimension(10) :: &
mma, &
mmb



!common/glradars/ndr,glrdr  !nochmal diskutieren
!common/ensemmtx/stks           !Elements of Backscatter and Extinction Matrices
!common/extmbmum/extm,bmum
include 'parameter.incf'
include 'constants.incf'

      idr = 1 ! initialize in case it gets skipped over later


       !LM gibt Temperatur in Kelvin aus, deshalb umrechnen
       ! Environment temperature in Celsius degrees


      temp=tmp-273.16

!      write(*,*) temp, 'temp_versuch'
      zhhx=0.
   thrdr(1)=0       !radar elevation angle, 0=horizontal
     nrank = 0

      !
      ! refractive index of air
      !

      refreair = 1
      refimair = 0


      !ndr: Number of Radar looking directions

      freq0=0.
      ndr = 1
      freq0=freq*1.e9
      wvln=vlght/freq0

!--------------------------------------------------
!set drop size distribution (DSD) for each category


  sznum=0.
  szlow=0.
  szupp=0.
  sdidv=0.
  sdgno=0.
  sdgmm=0.
  sdgml=0.
  lambd=0.

! Select the model hydrometeors types in the cloud, in the order of:
!     cldrp(cloud droplet),rains(rain drop),icclm(pristine ice column),
!     icplt(pristine ice plate),snows(snow flake),aggrs(aggregate),
!     graup(graupel),hlcnc(conical hail),hlsph(spheroidal hail).
!     e.g., hyddo='010000100' means to consider rain and graupel only.



  hyddo='000000000'
!  if (cloud > 0.) hyddo(1:1) = '1'
  if (rain > 0.)  hyddo(2:2) = '1'
  if (snow > 0.)  hyddo(4:4) = '1'
  if (graup > 0.) hyddo(7:7) = '1'

  tv=(temp+273.16)*(1. + 0.61*qd)       ! virtuelle Temperatur
  rhox=druck/(287.*tv)         ! Luftdichte über Gasgleichung
    call watereps(temp,epswtr,epsice)

!-->  cldrp
      if(hyddo(1:1) .eq. '1') then

         sznum(1)=8     !Number of quadrature points in DSD integration
         szlow(1)=0.001 !Lowest RAMS characteristic size (max or equi dia)
         szupp(1)=0.2   !Highest RAMS characteristic size (max or equi dia)
!         sdidv(1)=11    !Input device number for experimental DSD


         lambd(1) = (( rhox*cloud)/ (.82 * 1.E8 * (GAMMA(1.+3.)/GAMMA(1.))))**(1/(-1.-3))
         sdgno(1)=1000  !our N_0 in Gamma DSD
         sdgmm(1)=2     !RAMS mu=our m+1 in Gamma DSD
         sdgml(1)=.1/3.67 !RAMS D0=1/(our lambda) in Gamma DSD
      endif



!-->  rain
      if(hyddo(2:2) .eq. '1') then
         lambd(2) = 0.
         sznum(2)=13   !Number of quadrature points in DSD integration
         szlow(2)=0.01 !Daten aus DSD_Param *2 weil angabe in Radius

         sdgno(2) = 8.E6
!         sdgno(2) = 3.96e7 ! yen COSMO-DE

         sdgmm(2) = 0. ! mu-value in gamma dsd
!         sdgmm(2)=0.5  ! yen COSMO-DE

         lambd(2) = ((1000. * pi * sdgno(2))/(rhox * rain))**0.25    ! [m]
!         lambd(2) =  ((pi/6. * 1000. * sdgno(2) * GAMMA(4.+sdgmm(2))) / &
!            (rhox*rain))**(1./(4+sdgmm(2))) ! yen COSMO-DE

         sdgml(2) = 3.76/(lambd(2)*1.E-3)    ! mm
!         sdgml(2) = (3.67+sdgmm(2))/(lambd(2)*1.E-3)  ! yen COSMO-DE

         !  force upper limit of size integration to be 2.50*Do.
         !szupp(2) = 2.50 * sdgml(2)
         
         szupp(2) = min(2.5 * sdgml(2),5.0)


         sdidv(2)=12   !Input device number for experimental DSD data
!        The gamma parameters are input in fort.12 as a table and
!        read in dsdparam.f routine.
      endif


!-->  snow/icplt
      if(hyddo(4:4) .eq. '1') then
         lambd(4) = 0.
         sznum(4)=30    !Number of quadrature points in DSD integration
         szlow(4)=.001 !0.001Lowest RAMS characteristic size (max or equi dia)
         sdgmm(4)=0.   !2.     !RAMS mu=our m+1 in Gamma DSD
         
         
! yen 2008.12.09, field et al, dirty and quick
         !sdgno(4) = 8.E5
         !lambd(4) =  ((2*0.038*8.E5)/(rhox * snow))**(1./3.)
         !sdgml(4) = 3.67/(lambd(4)*1.E-3)    ! m Email von Th. Reinhardt
         !szupp(4) = 2.50 * sdgml(4)
         
! yen 2008.12.09, field et al, dirty and quick
  n0s1 = 13.5 * 5.65e5 ! parameter in N0S(T)
  n0s2 = -0.107        ! parameter in N0S(T), Field et al
  
  zams  = 0.069  ! Formfactor in the mass-size relation of snow particles corresponding to the radius

  mma = (/   5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
             0.312550,  0.000204,  0.003199, 0.000000, -0.015952 /)
  mmb = (/   0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
             0.060366,  0.000079,  0.000594, 0.000000, -0.003577 /)

  ! Calculate n0s using the tmerature-dependent moment
  ! relations of Field et al. (2005)
  tm = MAX(MIN(temp-273.15,0.0),-40.0)

  nn  = 3
  hlp = mma(1)      +mma(2)*tm      +mma(3)*nn       +mma(4)*tm*nn+mma(5)*tm**2 &
      + mma(6)*nn**2+mma(7)*tm**2*nn+mma(8)*tm*nn**2+mma(9)*tm**3+mma(10)*nn**3
  alf = 10.0**hlp
  bet = mmb(1)      +mmb(2)*tm      +mmb(3)*nn       +mmb(4)*tm*nn+mmb(5)*tm**2 &
      + mmb(6)*nn**2+mmb(7)*tm**2*nn+mmb(8)*tm*nn**2+mmb(9)*tm**3+mmb(10)*nn**3
  m2s = rhox*snow/zams
  m3s = alf*EXP(bet*LOG(m2s))

  hlp  = n0s1*EXP(n0s2*tm)
  n0s = 13.50 * m2s**4 / m3s**3
  n0s = MAX(n0s,0.5*hlp)
  n0s = MIN(n0s,1e2*hlp)
  n0s = MIN(n0s,1e9)
  n0s = MAX(n0s,1e6)

  lambda = ((2.*zams*n0s)/(rhox*snow))**(1./3.)

  sdgno(4) = n0s
  lambd(4) = lambda
  sdgml(4) = 2./(lambd(4)*1.e-3)
  szupp(4) = 4.*sdgml(4)
! yen 2008.12.09, field et al, dirty and quick

         

         thetas = 20.
         thetasd =10.
         epsair=cmplx(refreair,refimair)**2
!        density should be in  g/cc
     endif

!-->  graup
      if(hyddo(7:7) .eq. '1') then
         lambd(7) = 0.
         sdgno(7) = 0.
         sdgno(7) = 4.E6

! yen 2008.12.09 after mario mwmod, dirty and quick
         !lambd(7) =   ((pi*200.*sdgno(7))/(rhox*graup))**0.25
         !sdgml(7) = 3.76/(lambd(7)*1.E-3)    ! m
         !szupp(7) = 2.50 * sdgml(7)
         !if (szupp(7) .gt. 36) szupp(7) = 36
! yen 2008.12.09 after mario mwmod, dirty and quick
  lambd(7) = ((3.1*2.1*0.9135*169.9*sdgno(7))/(rhox*graup))**(1./4.1)         
  sdgml(7) = 3.1/lambd(7)
  szupp(7) = 4.*sdgml(4)
  if (szupp(7) .gt. 36) szupp(7) = 36 
! yen 2008.12.09 after mario mwmod, dirty and quick

         sznum(7)=10     !Number of quadrature points in DSD integration
         szlow(7)=0.1    !Lowest RAMS characteristic size (max or equi dia)
         sdidv(7)=17     !Input device number for experimental DSD data
         sdgmm(7)= 0. !1    !RAMS mu=our m+1 in Gamma DSD
!        density should be g/cc
         gradens=.2
!        compute the refractive index of graupel using Maxwell-Garnet formula
!        the parameter qg from RAMS will give the fraction of water in the
!        graupel
     qg=0
         thetag = 40
         thetagd = 20
         epsair=cmplx(refreair,refimair)**2
         call watereps(temp,epswtr,epsice)
         scmixgr=gradens/0.9
         if (temp .gt. 0.) then
           qg = 0.36
           !qg = 0.2  ! yen 2008.09.29
            thetag = 60
            thetagd = 40
          else
          qg = 0.
       endif
         if(qg .le. 0) then
            call refeffect(refregr,refimgr,scmixgr,epsair,epsice)
         elseif((qg .gt. 0) .and. (qg .le. .8))then
            scmixwater=qg/.8
            if(scmixwater .ge. (1-scmixgr))then
!              this condition implies that scmixair is zero
!              and set scmixwater is 1-scmixgr
               scmixwater = 1-scmixgr
               call refeffect(refregr,refimgr,scmixwater,epsice,epswtr)
            else
!              at this point scmixair is 1-scmixgr-scmixwater
               scmixwatermod=scmixwater/(scmixwater+scmixgr)
               call refeffect(refrep,refimp,scmixwatermod,epsice,epswtr)
               epsp=cmplx(refrep,refimp)**2
!              this is the intermediate effective dielectric constant of the
!              water/ice mixture
!              now calculate effective dielectric constant of air/water/ice mix
               call refeffect(refregr,refimgr,(scmixgr+scmixwater),epsair,epsp)
               !epsp ist eps von graupel, also ohne luft
            endif
         elseif(qg .gt. .8) then
!           at this point the graupel is all water
            refregr=real(csqrt(epswtr))
            refimgr=aimag(csqrt(epswtr))
         endif
      endif

!-->  hlcnc
      if(hyddo(8:8) .eq. '1') then
!         sdtyp(8)=2   !DSD type, 1=experimental, 2=Gamma.
         sznum(8)=25  !Number of quadrature points in DSD integration
         szlow(8)=0.1 !Lowest RAMS characteristic size (max or equi dia)
         szupp(8)=50. !Highest RAMS characteristic size (max or equi dia)
         sdidv(8)=18  !Input device number for experimental DSD data
         sdgno(8)=10  !our N_0 in Gamma DSD
         sdgml(8)=2   !RAMS D0=1/(our lambda) in Gamma DSD
         sdgmm(8)=2   !RAMS mu=our m+1 in Gamma DSD

!        density should be g/cc
         haildens=.89
     qh=.5
!         temphail=sctmp
         call watereps(temp,epswtr,epsice)
!        compute the refractive index of hail using Maxwell-Garnet formula
         if(qh .le. 0.)then
            refreha=real(csqrt(epsice))
            refimha=aimag(csqrt(epsice))
         endif
         if((qh .gt. 0.) .and. (qh .le. .8))then
            scmixhawater=qh/.8
            scmixhaice=1-scmixhawater

            call refeffect(refrewi,refimwi,scmixhaice,epswtr,epsice)
            epswi=cmplx(refrewi,refimwi)**2

            call refeffect(refreiw,refimiw,scmixhawater,epsice,epswtr)
            epsiw=cmplx(refreiw,refimiw)**2

            titu=0.2
            qitu=0.25

            xx= (((1-scmixhawater)/scmixhawater)-titu)/2*qitu
 !           xx=(((1-scmixhawater)/scmixhawater)-titu)/2*q
            erfct=erf(xx)

            epsreha=.5*((1-erfct)*real(epswi)+(1+erfct)*real(epsiw))
            epsimha=.5*((1-erfct)*aimag(epswi)+(1+erfct)*aimag(epsiw))
            epshail=cmplx(epsreha,epsimha)
            refreha=real(csqrt(epshail))
            refimha=aimag(csqrt(epshail))
         endif
         if(qh .gt. .8)then
            refreha=real(csqrt(epswtr))
            refimha=aimag(csqrt(epswtr))
         endif

      endif


!-->  hlsph
      if(hyddo(9:9) .eq. '1') then
!         sdtyp(9)=2   !DSD type, 1=experimental, 2=Gamma.
         sznum(9)=10  !Number of quadrature point in DSD integration
         szlow(9)=.1  !Lowest RAMS characteristic size (max or equi dia)
     szupp(9)=5  !Highest RAMS characteristice size (max or equi dia)
         sdidv(9)=19  !Input device number for experimental DSD data
         sdgno(9)=1735000  !our N_0 in Gamma DSD
         sdgml(9)=0.7008  !RAMS D0=1/(our lambda) in Gamma DSD
         sdgmm(9)=1   !RAMS mu=our m+1 in Gamma DSD

!        density should be g/cc
         haildens=.05
     qh=0.0
         epsair=cmplx(refreair,refimair)**2
!         tmphail=sctmp
!        the refractive index of sphericl hail is the same as conical hail
         call watereps(temp,epswtr,epsice)

         scmixha=haildens/0.9


         if(qh .le. 0)then

            call refeffect(refreha,refimha,scmixha,epsair,epsice)

!         if(qh .le. 0.)then
!            refreha=real(csqrt(epsice))
!            refimha=aimag(csqrt(epsice))
         endif
         if((qh .gt. 0.) .and. (qh .le. 0.80))then

            scmixhawater=qh/.8
            scmixhaice=1-scmixhawater
            call refeffect(refrewi,refimwi,scmixhaice,epswtr,epsice)
            epswi=cmplx(refrewi,refimwi)**2
            call refeffect(refreiw,refimiw,scmixhawater,epsice,epswtr)
            epsiw=cmplx(refreiw,refimiw)**2
            titu=0.2
            qitu=0.25
            xx=(((1-scmixhawater)/scmixhawater)-titu)/2*qitu
            erfct=erf(xx)
            epsreha=.5*((1-erfct)*real(epswi)+(1+erfct)*real(epsiw))
            epsimha=.5*((1-erfct)*aimag(epswi)+(1+erfct)*aimag(epsiw))
            epshail=cmplx(epsreha,epsimha)
            refreha=real(csqrt(epshail))
            refimha=aimag(csqrt(epshail))
         endif
         if(qh .gt. 0.80)then

            refreha=real(csqrt(epswtr))
            refimha=aimag(csqrt(epswtr))
         endif
   endif


extm = 0.
bmum = 0.


!
!===> Loop through the 9 hydrometeor catagories provided by RAMS model
!     to compute radar observables at a fixed space point.

      do 80 ist=1,9    !ist refers to our model hydrometeor type.
      if(hyddo(ist:ist).ne.'y' .and. hyddo(ist:ist).ne.'1') goto 80

      n0rank=5
      n0pnt=32
      crtnm=1.e-10

      arg1=sznum(ist)
      arg2=7*abs(sqrt(epswtr))*szupp(ist)*pi/wvln

      sznum(ist)=min(arg1,arg2)
      sznum(ist)=amax1(sznum(ist),8.0)


       nsz=INT(sznum(ist)+0.1)

      call gauslegquads(sz,wt,szlow(ist),szupp(ist),msz,nsz)




      do 70 iisz=1,nsz

    sdprb(iisz) = sdgno(ist)* (sz(iisz) *1.E-3) **sdgmm(ist)*exp(-lambd(ist)*sz(iisz)*1.E-3)

     dsz=sz(iisz)
      sdprb(iisz) = sdprb(iisz) * 1.E-3

       if (ist.eq.4) then
       pltdens = (6. * 0.038)/ (pi * 0.8 *sz(iisz))
       if (pltdens .gt. .2) pltdens = 0.2
       scmixplt=pltdens/0.9
       if (temp .gt. 0. ) then
          qs = 0.36
          !qs = 0.2   ! yen 2008.09.29
          thetas = 60.
          thetasd = 45.
          else
          qs = 0.
       endif

       if (qs .le. 0.) then
       call refeffect(refreplt,refimplt,scmixplt,epsair,epsice)
       endif


    if((qs .gt. 0.) .and. (qs .le. .8))then
            scmixwater=qs/.8
               if(scmixwater .ge. (1-scmixgr))then
!               write(*,*) 'scmixwater, scmixgr'
!               write(*,*) scmixwater, scmixgr
!              this condition implies that scmixair is zero
!              and set scmixwater is 1-scmixgr

            scmixwater=1-scmixplt
            call refeffect(refreplt,refimplt,scmixwater,epsice,epswtr)
            else
            scmixwatermod=scmixwater/(scmixwater+scmixplt)

            call refeffect(refrep,refimp,scmixwatermod,epsice,epswtr)
            epsp=cmplx(refrep,refimp)**2
            call refeffect(refreplt,refimplt,(scmixplt+scmixwater),epsair,epsp)
         endif
       endif
endif

      scmix=1.0


      if(ist.eq.1) call hydcldrp(dsz,epswtr,aovrb)
      if(ist.eq.2) call hydrains(dsz,epswtr,aovrb) !bestimmung der cadparameter ok
      if(ist.eq.4) call hydicplt(dsz,epsice,refreplt,refimplt, aovrb, thetas, thetasd)
      if(ist.eq.7) call hydgraup(dsz,epswtr,epsice,scmixgr,refregr,refimgr,aovrb, thetag, thetagd)
      if(ist.eq.8) call hydhlcnc(dsz,epswtr,epsice,scmix,refreha,refimha, aovrb)
      if(ist.eq.9) call hydhlsph(dsz,epswtr,epsice,scmix,refreha,refimha,aovrb)
!      write(*,*) refre, refim

!===> Generate canting parameters for integration over the normalized
!     particle orientation distribution.  This is processed for each
!     particle to reduce memory requirement.  These parameters will be
!     re-generate only if canting habit is changes from the previous one
!     as turned on by non-zero cdtyp parameter.
!

      call cadparam(iisz, thrdr)

!
!
!===> Prepare constant values defined in boby FSA, exp(+jwt) convention,
!     for computing scatterering amplitudes at any incidence directions
!     determined by particle canting to transform from boby frame to
!     radar frame.  e.g., Generate T-matrices of a spheroidal scatterer.
!

      iscshp=INT(scshp+0.01)

      arg1=dmx
      arg2=deq*aovrb**(2./3)/MIN(1.d0,aovrb)

      szpar=refre*max(arg1,arg2)/wvln/1000
!      write(*,*) ist
!     write(*,*) iscshp, aovrb, szpar, dmx, deq
!
!-->  use Mie solution for spherical particles, such as cloud droplet.
!
!***************************************************************
!     changed threshold for Mie, axis ratio must be .995 to 1.005
!     for spheres
!     aovrb verändert sich und dadurch nicht mehr miestreuung relevant

      if(iscshp.eq.1 .and. abs(aovrb-1.0).lt.0.005) then
!**************************************************************

        call miesphere()
        sctmodel='mie'
!        write(*,*) 'mie'
        goto 50
      endif


!
!-->  use Rayleigh Approximation for small particles, e.g., ice crystals
!
!  Next if statement commented out so that Rayleigh solution is used if
!   Mie is not.
     if(iscshp.eq.1 .and. szpar.lt.0.10) then

        call rayleighsph(aovrb)
!        write(*,*) 'rayleigh'
        sctmodel='ray'
        goto 40
      endif

!
!-->  use T-matrix method for any others.
!      write(*,*) 'tmat'

!     write(*,*) n0rank,n0pnt,crtnm,aovrb
      call tmtdriver(n0rank,n0pnt,crtnm,aovrb, nrank)
      sctmodel='tmt'
!
!===> Integration over the normalized particle orientation distribution.
!

!     write(*,*) 'aha'
40    call cadintegr(sctmodel, nrank)
!    stop
!     write(*,*) stks
50    idr=1
      extm(1,1,idr)=extm(1,1,idr)+sdprb(iisz)*stks( 1,idr)
      extm(2,1,idr)=extm(2,1,idr)+sdprb(iisz)*stks( 2,idr)
      extm(3,1,idr)=extm(3,1,idr)+sdprb(iisz)*stks( 3,idr)
      extm(4,2,idr)=extm(4,2,idr)+sdprb(iisz)*stks( 4,idr)
      extm(4,3,idr)=extm(4,3,idr)+sdprb(iisz)*stks( 5,idr)

      bmum(1,1,idr)=bmum(1,1,idr)+sdprb(iisz)*stks( 6,idr)
      bmum(2,1,idr)=bmum(2,1,idr)+sdprb(iisz)*stks( 7,idr)
      bmum(3,1,idr)=bmum(3,1,idr)+sdprb(iisz)*stks( 8,idr)
      bmum(4,1,idr)=bmum(4,1,idr)+sdprb(iisz)*stks( 9,idr)
      bmum(2,2,idr)=bmum(2,2,idr)+sdprb(iisz)*stks(10,idr)
      bmum(3,2,idr)=bmum(3,2,idr)+sdprb(iisz)*stks(11,idr)
      bmum(3,3,idr)=bmum(3,3,idr)+sdprb(iisz)*stks(12,idr)
      bmum(4,3,idr)=bmum(4,3,idr)+sdprb(iisz)*stks(13,idr)
      bmum(4,4,idr)=bmum(4,4,idr)+sdprb(iisz)*stks(14,idr)


70    continue
80    continue

!
!===> Fill up the uncomputed elements using the symmetry properties of
!     extinction and backscattering Mueller matrices.
!

      extm(1,1,idr)= extm(1,1,idr)*wvln
      extm(2,1,idr)= extm(2,1,idr)*wvln
      extm(3,1,idr)= extm(3,1,idr)*wvln
      extm(4,2,idr)= extm(4,2,idr)*wvln
      extm(4,3,idr)= extm(4,3,idr)*wvln
      extm(4,1,idr)= 0
      extm(1,2,idr)= extm(2,1,idr)
      extm(2,2,idr)= extm(1,1,idr)
      extm(3,2,idr)= 0
      extm(1,3,idr)= extm(3,1,idr)
      extm(2,3,idr)= 0
      extm(3,3,idr)= extm(1,1,idr)
      extm(1,4,idr)= 0
      extm(2,4,idr)=-extm(4,2,idr)
      extm(3,4,idr)=-extm(4,3,idr)
      extm(4,4,idr)= extm(1,1,idr)

      bmum(1,2,idr)= bmum(2,1,idr)
      bmum(4,2,idr)= 0
      bmum(1,3,idr)= bmum(3,1,idr)
      bmum(2,3,idr)= bmum(3,2,idr)
      bmum(1,4,idr)=-bmum(4,1,idr)
      bmum(2,4,idr)= 0
      bmum(3,4,idr)=-bmum(4,3,idr)


!
!===> Construct radar observables.
!
      zcoef=2.0e18/((pi/wvln)**4*abs((epswtr-1)/(epswtr+2))**2)

      call rdrobservable(zcoef,zhhx,ahhx,zvvx, zhvx,zvhx,delahvx, kkdp, rhohvx)




End subroutine versuch


    subroutine rdrobservable(zcoef,zhhx,ahhx,zvvx, zhvx, zvhx,delahvx, kkdp, rhohvx)



    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'


!      common/extmbmum/extm(4,4,mdr),bmum(4,4,mdr)
      common/rainrate/rain_gam,water_gam

!    common/glrdrobs/zhh(mdr),zvv(mdr),zhv(mdr),zvh(mdr),zdr(mdr),&
!        ldrh(mdr),ldrv(mdr),dph(mdr),dpv(mdr),rhohv(mdr),delhv(mdr),&
!        cdrl(mdr),orttl(mdr),aldl(mdr),dpl(mdr),cdrr(mdr),orttr(mdr),&
!        aldr(mdr),dpr(mdr),zhhdb(mdr),&
!        kdp(mdr),ahh(mdr),avv(mdr),delahv(mdr)

      real dph1(mdr), dph2(mdr), dph3(mdr),dph4(mdr), dpv1(mdr), &
    dpv2(mdr), dpv3(mdr), dpv4(mdr), rhohv1(mdr), rhohv2(mdr),&
    rhohv3(mdr), rhohv4(mdr), w11, w21, w22, w12,w31r,w31i,w32r,&
        w32i, yy, xx, zhhx,ahhx, zzzdr,llldr, delahvx, kkdp, &
                rhohvx, zvvx, zhvx, zvhx
    real*8 :: degrad
    integer ::  numcall, idr




    include 'constants.incf'

      data degrad /1.745329251994329e-2/

!      data numcall /0/        ! Number of calls


!===> Calculate polarimetric radar parameters
!
      numcall = 0
      idr=1
!
!===> Backscattering related radar parameters:
!
!-->  linear polarization:
!     zhh,zvv,zvh,zhv:effective reflectivity factor, mm^6/m^3
!     zdr: differentail reflectivity, dB
!     ldrh,ldrv: linear depolarization ratio, dB
!     dph,dpv: linear degree of polarization, %
!     rhohv: zero-copolar cross-correlation coefficient, no unit
!     delhv: backscattering differential phase shift, degree
!
!         write(*,*) bmum(1,1,idr)
!stop
      zvv(idr)=zcoef*(bmum(1,1,idr)+bmum(2,1,idr)+bmum(1,2,idr)+bmum(2,2,idr))
!      write(*,*) 'zvv', zvv
!stop
!      write(*,*) zvv(idr), zcoef, bmum(1,1,idr)
      zvh(idr)=zcoef*(bmum(1,1,idr)+bmum(2,1,idr)-bmum(1,2,idr)-bmum(2,2,idr))
      zhv(idr)=zcoef*(bmum(1,1,idr)-bmum(2,1,idr)+bmum(1,2,idr)-bmum(2,2,idr))
!    write(*,*) bmum(1,1,idr)
      zhh(idr)=zcoef*(bmum(1,1,idr)-bmum(2,1,idr)-bmum(1,2,idr)+bmum(2,2,idr))

      zhhdb(idr)=fdb(zhh(idr))

      zdr(idr)=0.

      zdr(idr)=fdb(zhh(idr)/zvv(idr))


     ! if (zhh(idr) .lt. 1.E-6) then
     !    zdr(idr) = 0.
     ! endif

     ! if (zvv(idr) .lt. 1.E-6) then
     !    zvv(idr) = 0.
     ! endif
!      write(*,*) 'zhv, zhh'
!      write(*,*) zhv(idr), zhh(idr)
      ldrv(idr)=fdb(zvh(idr)/zvv(idr))
      ldrh(idr)=fdb(zhv(idr)/zhh(idr))




!neu
    dph1(idr)=(bmum(2,1,idr)-bmum(2,2,idr))**2
    dph2(idr)=(bmum(3,1,idr)-bmum(3,2,idr))**2
    dph3(idr)=(bmum(4,1,idr)-bmum(4,2,idr))**2
    dph4(idr)= bmum(1,1,idr)-bmum(1,2,idr)
        dph(idr)=100*sqrt(dph1(idr)+dph2(idr)+dph3(idr))/dph4(idr)



    dpv1(idr)= (bmum(2,1,idr)+bmum(2,2,idr))**2
    dpv2(idr)= (bmum(3,1,idr)+bmum(3,2,idr))**2
    dpv3(idr)= (bmum(4,1,idr)+bmum(4,2,idr))**2
    dpv4(idr)= bmum(1,1,idr)+bmum(1,2,idr)
        dpv(idr)=100*sqrt(dpv1(idr)+dpv2(idr)+dpv3(idr))/dpv4(idr)



    rhohv1(idr)=(bmum(3,3,idr)+bmum(4,4,idr))**2
    rhohv2(idr)=(bmum(4,3,idr)-bmum(3,4,idr))**2
    rhohv3(idr)= bmum(1,1,idr)-bmum(2,1,idr)-bmum(1,2,idr)+bmum(2,2,idr)
    rhohv4(idr)= bmum(1,1,idr)+bmum(2,1,idr)+bmum(1,2,idr)+bmum(2,2,idr)
        rhohv(idr)=sqrt((rhohv1(idr)+rhohv2(idr))/(rhohv3(idr)*rhohv4(idr)))




      delhv(idr)=1.0/degrad*atan2((bmum(4,3,idr)-bmum(3,4,idr)),(bmum(3,3,idr)+bmum(4,4,idr)) )

!
!-->  circular polarization:
!     cdrl,cdrr: circular depolarization ratio, dB
!     corttl,orttr: apparent degree of orientation, no unit
!     caldl,aldr: apparent orientation angle, degree
!     cdpl,dpr: circular degree of polarization, %

      w11=0.5*(bmum(1,1,idr)+bmum(1,4,idr)-bmum(4,1,idr)-bmum(4,4,idr))
      w21=0.5*(bmum(1,1,idr)+bmum(1,4,idr)+bmum(4,1,idr)+bmum(4,4,idr))
      w22=0.5*(bmum(1,1,idr)-bmum(1,4,idr)+bmum(4,1,idr)-bmum(4,4,idr))
      w12=0.5*(bmum(1,1,idr)-bmum(1,4,idr)-bmum(4,1,idr)+bmum(4,4,idr))
      w31r= 0.5*(bmum(2,1,idr)+bmum(2,4,idr))
      w31i=-0.5*(bmum(3,1,idr)+bmum(3,4,idr))
      w32r= 0.5*(bmum(2,1,idr)-bmum(2,4,idr))
      w32i= 0.5*(bmum(3,1,idr)-bmum(3,4,idr))
      cdrl(idr)=fdb(w11/w21)
      cdrr(idr)=fdb(w22/w12)
      orttl(idr)=sqrt((w31r*w31r+w31i*w31i)/(w11*w21))
      orttr(idr)=sqrt((w32r*w32r+w31i*w31i)/(w12*w22))
      aldl(idr)=0.5*(atan(w31i/w31r)-pi)/degrad
      aldr(idr)=0.5*(atan(w32i/w32r)-pi)/degrad

      dpl(idr)=100*sqrt((bmum(2,1,idr)+bmum(2,4,idr))**2+(bmum(3,1,idr)+bmum(3,4,idr))**2&
          +(bmum(4,1,idr)+bmum(4,4,idr))**2)/(bmum(1,1,idr)+bmum(1,4,idr))
      dpr(idr)=100*sqrt((bmum(2,1,idr)-bmum(2,4,idr))**2+(bmum(3,1,idr)-bmum(3,4,idr))**2&
         +(bmum(4,1,idr)-bmum(4,4,idr))**2)/(bmum(1,1,idr)-bmum(1,4,idr))

!
!==> Forward scattering related radar parameters:
!    ahh,avv: specific attenuation, db/km
!    kdp: specific differential phase, degree/km
!
      avv(idr)=8.686e3*0.5*(extm(1,1,idr)+extm(2,1,idr))
      ahh(idr)=8.686e3*0.5*(extm(1,1,idr)-extm(2,1,idr))
      delahv(idr)=ahh(idr)-avv(idr)
      kdp(idr)=extm(3,4,idr)/degrad*1000
      kkdp =  kdp(idr)

      yy=gammln(4.1)
      xx=erf(6.25e-3)
      zzzdr = zdr(1)


      llldr = ldrh(1)
      rhohvx = rhohv(1)
!      write(*,*) 'LDR', llldr
!*****************************************************************
!    changed order of variables as shown below
      numcall = numcall + 1
! WY 2007.08.30
      !if (numcall .eq. 1) zhhx=fdb(zhh(1))
      zhhx = zhh(1)
! WY 2007.08.30
      zhvx = zhv(1)
      zvvx = zvv(1)
      zvhx = zvh(1)
 !     write(*,*) zhhx
      ahhx = ahh(1)
      delahvx = delahv(1)
!write(55,*) '   zhh   zdr   kdp    ahh       delahv'&
!     &     //'    rhohv   delhv     LDR  rainrate'
!    If you want to output rainrate then change water_gam to rain_gam
!      write(55,44) fdb(zhh(1)),zdr(1),kdp(1),ahh(1),delahv(1),rhohv(1)&
!     ,delhv(1),ldrh(1),rain_gam


 44   format(f9.3,2f7.3,2(1x,E10.3),f8.4,3f8.2)
      return
      end subroutine







!=======================================================================
!
      subroutine cadparam(isz, thrdr)



!  use funktion



    implicit none
!
!-----------------------------------------------------------------------
!     This routine returns the canting parameters for use in the canting
!     angle distribution integration.
!
!-->  Orientation of a rotationally and equatorially symmetric particle
!     is sufficiently defined by two canting angles, the polar (1st) and
!     the azimuthal (2nd) angles of its symmetry or body Z-axis.
!
!     Array glcad stores the control parameters for each Canting Angle
!     Distribution in the order of:  CAD type, lower limit, upper limit,
!     number of quadrature points, mean value, standard deviation.
!     Array glrdr stores the polar and azimuthal angles of radar beams.
!
!     Currently supported CADs and parameters (in degrees) used are:
!     1=Uniform: cdtyp=1,calow,caupp,canum.
!     2=Gaussian: cdtyp=2,calow,caupp,canum,caavr,cadev.
!     3=Simple harmonic oscillator: cdtyp=3,calow,caupp,canum,caavr.
!
!-->  The following parameters can be fine tuned to speed up the code,
!
!     cdtyp (in hydxxxx routines): if any of the two cdtyp is turned off
!           (=<0), orientation habit at the previous step will be used,
!     canum (in hydxxxx routines): try to set them to the smallest,
!     cdpskp: skip current step if cdprb/cdpmx<cdpskp (this is extremely
!             helpful when canting distribution varies sharply),
!     cdpdth: approximate scattering amplitudes at incidence angles in
!             the range of [x-cdpdth/2, x+cdpdth/2) by that at angle x.
!-----------------------------------------------------------------------
!

      include 'variablen.incf'
      include 'parameter.incf'


!      common/glradars/ndr,glrdr(mdr,2)
!      common/glhydbdy/glbdy(12),glcad(2,6)
!

      integer ::  idr, icd, ica, ica1, ica2, icp, ndr

!      common/cantpara/cdpnum,cdprb,cdpth,&
!          cdpab,cdpdd,cdpab2,cdpdd2,cdpabd
!

      integer :: ncdp, jcdp, icdp, isz,n
      real, parameter  :: cdpskp=1.e-5,cdpdth=0.1
      integer, parameter :: mdt=1801   !mdt=180/cdpdth+1
      integer ncp(mdr)
      real  pth(mdt),prb(mdt),pab(mdt),pdd(mdt),pab2(mdt),pdd2(mdt),&
            pabd(mdt), runit(3),txr(3),teu(3,3),trdr(3,3,mdr),wt(mca),&
         cpprb(mca,mca),prbca(mca,2),cang(mca),sinca(mca,2),cosca(mca,2),&
         thcp(mca2,mdr),prcp(mca2,mdr),p1cp(mca2,mdr),p2cp(mca2,mdr)
      real thrdr(2)
      real, parameter :: degrad=1.745329251994329e-2
      real :: cdpdt0, cosph, cdnrm,sinth,costh,&
          sinph,cdpwt,ab,cc, eu2t,eu1t,&
          x,cdpmx,theui,tmp1,pheui,tmp2,tmp3,dd
      real ccanum(2),ccdtyp,ccalow, ccaavr, ccadev, ccaupp

      include 'constants.incf'
       ndr = 1


!
!===> Not to skip this routine if orientation habit has been changed
!
      if(isz.ne.1.and.(cdtyp(1).le.0.0)) goto 140
!      print*,isz,glcad(1,1),glcad(2,1)
!
!===> Calculate Transformation matrices from spherical to Cartesian
!     coordinates for all radar directions measured in the lab frame,
!
      cdpdt0=cdpdth*degrad

      do 10 idr=1,ndr
      sinth=sin((90-thrdr(idr))*degrad)
      costh=cos((90-thrdr(idr))*degrad)
      sinph=sin(glrdr(idr,2)*degrad)
      cosph=cos(glrdr(idr,2)*degrad)
      trdr(1,1,idr)= sinth*cosph
      trdr(2,1,idr)= sinth*sinph
      trdr(3,1,idr)= costh
      trdr(1,2,idr)= costh*cosph
      trdr(2,2,idr)= costh*sinph
      trdr(3,2,idr)=-sinth
      trdr(1,3,idr)=-sinph
      trdr(2,3,idr)= cosph
      trdr(3,3,idr)= 0.0d0
!     write(*,*) trdr ok!!
      ncp(idr)=0
 10   continue
!
!===> Generate (sin and cos) canting angles and normalized probabilities
!     at the nested Gaussian-Legendre quadrature grid points.
!     Note: Include the sin factor of the solid angle in the probability
!     of the first canting angle or the polar angle of the body z-axis.
!


      do 30 icd=1,2

      ccanum(icd)=canum(icd)+0.01
      ccdtyp=cdtyp(icd)+0.01
      ccalow=calow(icd)*degrad
      ccaupp=caupp(icd)*degrad
      ccaavr=caavr(icd)*degrad
      ccadev=cadev(icd)*degrad
!    write(*,*) caavr(1), cadev(1)
!stop
      if(ccanum(icd).eq.1) then
      cang(1)=0.5*(ccalow+ccaupp)
      sinca(1,icd)=sin(cang(1))
      cosca(1,icd)=cos(cang(1))
      prbca(1,icd)=1.0
      goto 30
      else
      n=ccanum(icd)
      call gauslegquads(cang,wt,ccalow,ccaupp,mca,n)
      do 20 ica=1,ccanum(icd)
!      print*,ica,ccanum(icd),'here'
      x=cang(ica)
      sinca(ica,icd)=sin(x)
      cosca(ica,icd)=cos(x)
      prbca(ica,icd)=wt(ica)*cadfunc(x,ccdtyp,ccalow,ccaupp,ccaavr,ccadev)
!      if(icd .eq. 1)print*,prbca(ica,icd),' here'
      if(icd.eq.1) prbca(ica,1)=prbca(ica,1)*sinca(ica,1)
 20   continue
      endif
 30   continue
!    write(*,*) canum ok

      cdnrm=0.0
      cdpmx=0.0
      do 40 ica2=1,ccanum(2)
      do 40 ica1=1,ccanum(1)
      cpprb(ica1,ica2)=prbca(ica1,1)*prbca(ica2,2)
      if(cpprb(ica1,ica2).gt.cdpmx) cdpmx=cpprb(ica1,ica2)
      cdnrm=cdnrm+cpprb(ica1,ica2)
 40   continue
      cdpmx=cdpmx/cdnrm
!
!===> Main loop over the two canting angles.  Skip current step if CAD
!     probability factor is less than cdpskp of its maximum.
!
      do 70 ica2=1,ccanum(2)
      teu(2,1)=-sinca(ica2,2)
      teu(2,2)= cosca(ica2,2)
      teu(2,3)= 0
      do 60 ica1=1,ccanum(1)
      cpprb(ica1,ica2)=cpprb(ica1,ica2)/cdnrm
      cdpwt=cpprb(ica1,ica2)
      if( (cdpwt/cdpmx).lt.cdpskp ) goto 60

      teu(1,3)=-sinca(ica1,1)
      teu(3,3)= cosca(ica1,1)
      teu(1,1)= teu(3,3)*teu(2,2)
      teu(3,1)=-teu(1,3)*teu(2,2)
      teu(1,2)=-teu(3,3)*teu(2,1)
      teu(3,2)= teu(1,3)*teu(2,1)
!
!===> Mapping the directions of interest from lab frame to body frame.
!     compute the polar and azimuthal angles measureed in body FSA and
!     the two transformation parameters from body FSA to lab FSA for the
!     transverse E-components at each incidence direction.
!
      do 50 idr=1,ndr

      runit(1)=teu(1,1)*trdr(1,1,idr)+teu(1,2)*trdr(2,1,idr)+teu(1,3)*trdr(3,1,idr)
      runit(2)=teu(2,1)*trdr(1,1,idr)+teu(2,2)*trdr(2,1,idr)
!    1        +teu(2,3)*trdr(3,1,idr)
      runit(3)=teu(3,1)*trdr(1,1,idr)+teu(3,2)*trdr(2,1,idr)+teu(3,3)*trdr(3,1,idr)

      theui=acos( runit(3) )
      tmp1=runit(1)*runit(1) + runit(2)*runit(2)
      if( tmp1.gt.1.e-14) then
      pheui=atan2( runit(2), runit(1) )
      else if( runit(3).gt.0.0 ) then
      pheui=0.0
      else
      pheui=pi
      endif

      txr(1)=cos(theui)
      txr(2)=txr(1)*sin(pheui)
      txr(1)=txr(1)*cos(pheui)
      txr(3)=-sin(theui)
!
!-->  compute the two transformation parameters from the body spherical
!     to the lab spherical coordinates for the transverse components of
!     each of the directions of interest.
!
      tmp1=txr(1)*teu(1,1)+txr(2)*teu(2,1)+txr(3)*teu(3,1)
      tmp2=txr(1)*teu(1,2)+txr(2)*teu(2,2)+txr(3)*teu(3,2)
      tmp3=txr(1)*teu(1,3)+txr(2)*teu(2,3)+txr(3)*teu(3,3)
      eu1t=trdr(1,2,idr)*tmp1+trdr(2,2,idr)*tmp2+trdr(3,2,idr)*tmp3
      eu2t=-( trdr(1,3,idr)*tmp1+trdr(2,3,idr)*tmp2 )

      ncp(idr)=ncp(idr)+1
      thcp(ncp(idr),idr)=theui
      prcp(ncp(idr),idr)=cdpwt
      p1cp(ncp(idr),idr)=eu1t
      p2cp(ncp(idr),idr)=eu2t

 50   continue

 60   continue
 70   continue
!
!===> Generate canting parameters
!
!-->  exact computation, as denoted by cdpdt0 =< 0
!
      if(cdpdt0.le.0.0) then
      do 90 idr=1,ndr

      cdpnum(idr)=ncp(idr)
!      write(*,*) cdpnum
!stop
      do 80 icp=1,ncp(idr)
      cdpth(icp,idr)=thcp(icp,idr)
      cc=prcp(icp,idr)
      ab=p1cp(icp,idr)*p2cp(icp,idr) *2
      dd=p1cp(icp,idr)*p1cp(icp,idr)-p2cp(icp,idr)*p2cp(icp,idr)
      cdprb(icp,idr)=cc
      cdpab(icp,idr)=cc * ab
      cdpdd(icp,idr)=cc * dd
      cdpab2(icp,idr)=cc * ab*ab
      cdpabd(icp,idr)=cc * ab*dd
      cdpdd2(icp,idr)=cc * dd*dd
 80   continue
 90   continue
      goto 140
      endif
!
!===> Approximate computation.  first sort incidence angle in intervals.
!
      ncdp=pi/cdpdt0+1

     do 130 idr=1,ndr

      do 100 icdp=1,ncdp
      prb(icdp)=0.0
      pab(icdp)=0.0
      pdd(icdp)=0.0
      pab2(icdp)=0.0
      pabd(icdp)=0.0
      pdd2(icdp)=0.0
 100  continue

      do 110 icp=1,ncp(idr)
!         write(*,*) icp, 'icp'
      icdp=thcp(icp,idr)/cdpdt0+0.5
      cc=prcp(icp,idr)
      ab=p1cp(icp,idr)*p2cp(icp,idr) *2
      dd=p1cp(icp,idr)*p1cp(icp,idr)-p2cp(icp,idr)*p2cp(icp,idr)
      pth(icdp)=(icdp-1)*cdpdt0
      prb(icdp)=prb(icdp) + cc
!      write(*,*) cc, 'cc'
      pab(icdp)=pab(icdp) + cc*ab
      pdd(icdp)=pdd(icdp) + cc*dd
      pab2(icdp)=pab2(icdp) + cc*ab*ab
      pabd(icdp)=pabd(icdp) + cc*ab*dd
      pdd2(icdp)=pdd2(icdp) + cc*dd*dd
 110  continue

      jcdp=0
      do 120 icdp=1,ncdp
!         write(*,*) icdp, prb(icdp), 'test3'
      if(prb(icdp).gt.0.0) then
      jcdp=jcdp+1
      cdpth(jcdp,idr)=pth(icdp)
      cdprb(jcdp,idr)=prb(icdp)
      cdpab(jcdp,idr)=pab(icdp)
      cdpdd(jcdp,idr)=pdd(icdp)
      cdpab2(jcdp,idr)=pab2(icdp)
      cdpabd(jcdp,idr)=pabd(icdp)
      cdpdd2(jcdp,idr)=pdd2(icdp)
      endif
!      print*,cdprb(jcdp,idr),cdpab(jcdp,idr),' here'
 120  continue
      cdpnum(idr)=jcdp
!
 130  continue
!
!
!     write(*,*) cdpnum(1), jcdp, 'test2'
!top
 140  return
!top
      end subroutine cadparam





       subroutine cadintegr(sctmodel, nrank)

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'


!
!-----------------------------------------------------------------------
!     This routine integrates the extinction and backscattering Mueller
!     matrice over the normalized particle orientation distributions.
!     Only independent elements are processed.
!-----------------------------------------------------------------------
!
      character sctmodel*3

!

      integer  ncdp,idr,icdp, nrank
!      common/cantpara/cdpnum(mdr),cdprb(mca2,mdr),cdpth(mca2,mdr),&
!         cdpab(mca2,mdr),cdpdd(mca2,mdr),cdpab2(mca2,mdr),&
!         cdpdd2(mca2,mdr),cdpabd(mca2,mdr)
!      common/ensemmtx/stks(14,mdr)
!
    real :: thinc, ab, ab2, dd1, dd2,cc,abd
      real sem(7)
    include 'constants.incf'
!
!===> Initialize array stks.
!
!        write(*,*) 'hurra'
       idr=1

!      do 10 i=1,14
      stks=0

! 10   continue

      ncdp=cdpnum(idr)
!      write(*,*) 'neu', ncdp, cdpnum
      do 40 icdp=1,ncdp
!write(*,*) icdp, 'test'
      thinc=cdpth(icdp,idr)
      cc=cdprb(icdp,idr)
      ab=cdpab(icdp,idr)
      dd1=cdpdd(icdp,idr)
      ab2=cdpab2(icdp,idr)
      dd2=cdpdd2(icdp,idr)
      abd=cdpabd(icdp,idr)
!    write(*,*) ab2, dd2, abd
!top
!stop
!      print*,cc,ab,'prob'
!
!===> Calculate scattering amplitudes in body FSA, exp(+jwt) convention,
!     and perform coordinate transformation, and calculate expectation.
!     Note: the order of array order bstks: svvsvv,shvshv,svhsvh,shhshh,
!     ReImsvvshv,ReImsvvsvh,ReImsvvshh,ReImshvsvh,ReImshvshh,ReImsvhshh;
!     and the order of array fstks is: ReImsvv,ReImshv,ReImsvh,ReImshh.
!


      call scatamp(thinc,sem,sctmodel, nrank)
!      write(*,*) sem
!stop
!
      stks( 1,idr)=stks( 1,idr) + cc*sem(1)

      stks( 2,idr)=stks( 2,idr) + dd1*sem(2)
      stks( 3,idr)=stks( 3,idr) - ab*sem(2)
      stks( 4,idr)=stks( 4,idr) + ab*sem(3)
      stks( 5,idr)=stks( 5,idr) + dd1*sem(3)

      stks( 6,idr)=stks( 6,idr) + cc*sem(4)
      stks( 7,idr)=stks( 7,idr) + dd1*sem(5)
      stks( 8,idr)=stks( 8,idr) +abd*(sem(4)+sem(6))
      stks( 9,idr)=stks( 9,idr) + ab*sem(7)
      stks(10,idr)=stks(10,idr) +dd2*sem(4)-ab2*sem(6)
      stks(11,idr)=stks(11,idr) + ab*sem(5)
      stks(12,idr)=stks(12,idr) + cc*sem(6)
      stks(13,idr)=stks(13,idr) + dd1*sem(7)
      stks(14,idr)=stks(14,idr) +dd2*sem(6)-ab2*sem(4)

 40   continue
!      print*,stks(1,idr),stks(2,idr),stks(3,idr),stks(4,idr),'stksi'

!
! 50   continue
!
!
!      write(*,*) stks
!stop
      return
      end subroutine


!=======================================================================
end module src_versuch
