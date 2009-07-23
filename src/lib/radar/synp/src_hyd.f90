module src_hyd
      
      use src_hilfe
      use src_global
      contains


              subroutine hydaggrs(dsz,epsice,scmix0,aovrb,refreagg,refimagg)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     aggregate of characteristic size dsz=deq
!-----------------------------------------------------------------------
!
      complex epsice
      real  dsz,glbdy(12),glcad(2,6), &
    tempenv,temprain,tempgraup,temphail,tempcld,&
                      tempclm,tempplt,tempsnow,tempagg, &
     scmixg,scmixh
!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!                      tempclm,tempplt,tempsnow,tempagg
      common/heatbudget/scmixg,scmixh
      real refreclm,refimclm,refresnow,refimsnow,refreagg, &
                refimagg,refregr,refimgr,refreha,refimha, &
                refreair,refimair,refreplt,refimplt
      
      real snowdens,aggdens,gradens,haildens,clmdens,pltdens
      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!         sctmp,scmix, scmix0
      real*8 aovrb
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!        (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!        (glbdy(7),deq),(glbdy(8),dmx), &
!        (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      refre=refreagg
      refim=refimagg
      scshp= 1    ! 1 denotes sphere/spheroid shape
      deq=dsz
      aovrb=1.0
      dmx=deq*aovrb**(2.d0/3.d0)/min(1.0d0,aovrb)
!
!===> Orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt. 1.00) cdtyp(1)=0
      canum(1)= 10
      calow(1)=  0
      caupp(1)= 90
      caavr(1)=  0
      cadev(1)= 45

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360


      return
      end subroutine

        subroutine hydgraup(dsz,epswtr,epsice,scmix0, refregr, refimgr, aovrb, thetag, thetagd)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     conical graupel of characterictic size dsz=wcnzc
!-----------------------------------------------------------------------
!
      complex epswtr,epsice
      real dsz,glbdy(12),glcad(2,6), &
!    tempenv,temprain,tempgraup,temphail,tempcld,&     
!          tempclm,tempplt,tempsnow,tempagg, &
      scmixg,scmixh
!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,& 
!            tempclm,tempplt,tempsnow,tempagg
!     common/heatbudget/scmixg,scmixh
      
      real refregr, refimgr
      real*8 aovrb
!      real snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!         sctmp,scmix, scmix0
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!        (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!        (glbdy(7),deq),(glbdy(8),dmx), &
!        (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      refre=refregr
      refim=refimgr
      scshp=1    !2 denotes Wang's Conical shape
      wcnzc=dsz
      wcnxl=dsz*0.9
      wcnlm=2
      aovrb=0.5 !0.5
      dmx=dsz
      deq=dmx/aovrb**(2.d0/3.d0)*min(1.0d0,aovrb)
!      deq = deq*((1./0.9)**(1./3.))
!
!===> Orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt.1.00) cdtyp(1)=0
      canum(1)= 20
      calow(1)=  0
      caupp(1)= thetag
      caavr(1)=  0
      cadev(1)= thetagd

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360


      return
      end subroutine

      subroutine hydhlcnc(dsz,epswtr,epsice,scmix0,refreha,refimha, aovrb)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     conical hail of characteristic size dsz=wcnxl.
!-----------------------------------------------------------------------
!
      complex epswtr,epsice
      real dsz, glbdy(12),glcad(2,6), &
    tempenv,temprain,tempgraup,temphail,tempcld,&
          tempclm,tempplt,tempsnow,tempagg, &
        scmixg,scmixh, scmix0
!      common/glhydbdy/glbdy(12),glcad(2,6)
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!          tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
      real refreclm,refimclm,refresnow,refimsnow,refreagg, &
            refimagg,refregr,refimgr,refreha,refimha, &
            refreair,refimair,refreplt,refimplt
       real*8 aovrb
!      real  snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!        sctmp,scmix
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!        (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!        (glbdy(7),deq),(glbdy(8),dmx), &
!        (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      refre=refreha
      refim=refimha
      scshp= 2    ! 2 denotes Wang's conical shape
      wcnxl=dsz
      wcnzc=dsz*0.8
      wcnlm=20
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt.1.00) cdtyp(1)=0
      canum(1)= 20
      calow(1)=  0
      caupp(1)=180
      caavr(1)=  0
      cadev(1)= 20

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360


      return
      end subroutine


            subroutine hydhlsph(dsz,epswtr,epsice,scmix0,refreha, refimha,aovrb)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     oblate hail of characterictic size dsz=deq.
!-----------------------------------------------------------------------
!
      complex epswtr,epsice
      real dsz, scmix0,glbdy(12),glcad(2,6) 
!      common/glhydbdy/glbdy,glcad

    real  tempenv,temprain,tempgraup,temphail,tempcld,&
        tempclm,tempplt,tempsnow,tempagg, &
    scmixg,scmixh
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!        tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
      real refreclm,refimclm,refresnow,refimsnow,refreagg, &
          refimagg,refregr,refimgr,refreha,refimha,  &
          refreair,refimair,refreplt,refimplt
      real*8 aovrb
!      real snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!       sctmp,scmix
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!        (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!        (glbdy(7),deq),(glbdy(8),dmx), &
!        (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
 !     equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      refre=refreha
      refim=refimha
      scshp= 1    ! 1 denotes sphere/spheroid shape
      deq=dsz

!      aovrb=1.00
      aovrb=0.8
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      canum(1)= 10
      calow(1)= 0
      caupp(1)= 90
      caavr(1)= 0
      cadev(1)= 45

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)= 0 
      caupp(2)= 360


      return
      end subroutine


      subroutine hydicclm(dsz,epsice,aovrb,refreclm,refimclm)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     colume-like pristine crystal of characterictic size dsz=dmx
!-----------------------------------------------------------------------
!
      complex epsice
      real glbdy(12),glcad(2,6), &
         tempenv,temprain,tempgraup,temphail,tempcld,&
         tempclm,tempplt,tempsnow,tempagg, &
         scmixg,scmixh, dsz
!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!                        tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
    real*4 refreclm,refimclm
      complex refresnow,refimsnow,refreagg, &
             refimagg,refregr,refimgr,refreha,refimha, &
             refreair,refimair,refreplt,refimplt
      common/refractive/refresnow,refimsnow,refreagg, &
             refimagg,refregr,refimgr,refreha,refimha, &
             refreair,refimair,refreplt,refimplt
!      real  snowdens,aggdens,gradens,haildens,clmdens,pltdens, &
!    real    refre, refim
      real*8 aovrb, rrefre, rrefim
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!         sctmp,scmix,scshp,deq,dmx,wcnxl,wcnzc,wcnlm

!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!       (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!       (glbdy(7),deq),(glbdy(8),dmx), &
!       (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      sctmp=tempenv
!      call watereps(freq,sctmp,epswtr,epsice)
!      refre= real(sqrt(epsice))
!      refim=aimag(sqrt(epsice))
    refre = refreclm
    refim = refimclm
      

      scshp= 1     !1 denotes sphere/spheroid shape
      dmx=dsz
      aovrb = 1.0
!      aovrb=10
      deq=dmx/aovrb**(2.d0/3.d0)*min(1.0d0,aovrb)
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt.0.07) cdtyp(1)=0
      canum(1)= 10
      calow(1)=  0
      caupp(1)= 90
      caavr(1)= 90
      cadev(1)= 20

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360


      return
      end subroutine


      subroutine hydicplt(dsz,epsice,refreplt, refimplt, aovrb, thetas, thetasd)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     plate-like pristine crystal of characteristic size dsz=dmx
!-----------------------------------------------------------------------
!
      complex epsice, epswtr
      real  dsz,glbdy(12),glcad(2,6) , &
        tempenv,temprain,tempgraup,temphail,tempcld,&
           tempclm,tempplt,tempsnow,tempagg, &
       scmixg,scmixh

!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!           tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
      real refreclm,refimclm,refresnow,refimsnow,refreagg, &
           refimagg,refregr,refimgr,refreha,refimha, &
           refreair,refimair,refreplt,refimplt, thetas, thetasd
      real*8 aovrb
!      real snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2),  freq, &
!           sctmp,scmix, hoehe
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!       (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!       (glbdy(7),deq),(glbdy(8),dmx), &
!       (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
 !     (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!      (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      refre=refreplt
      refim=refimplt
      
      
      scshp=1    ! 1 denotes sphere/spheroid shape
      dmx=dsz
      aovrb =0.3
!      hoehe = 1.3E-2 * (2* dsz)**0.375
!      hoehe = 0.138 * dsz**0.778     !Pruppacher S. 51, ohne 0.5
!      aovrb = (hoehe/dsz)
!      aovrb = 0.1
!      aovrb=0.1    ! should size dependent, consult with Bob Walco

      deq=dmx/aovrb**(2.d0/3.d0)*min(1.0d0,aovrb)
     !  deq = dmx !deq equivolume diameter, dmx maximum dimension
!      write(*,*) dmx, deq 
!      deq = deq * ((1./0.6)**(1./3.)) !volumenkorrektur günther
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt. 1.00) cdtyp(1)=0
      canum(1)= 10
      calow(1)=  0
      caupp(1)= thetas
      caavr(1)=  0
      cadev(1)= thetasd

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360


      return
      end subroutine

 subroutine hydcloudice(dsz,epsice,refreclice,refimclice,aovrb)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     snow flake of characteristic size dsz=dmx
!-----------------------------------------------------------------------
!
      complex epsice
      real  dsz, glbdy(12),glcad(2,6), &
    tempenv,temprain,tempgraup,temphail,tempcld,&
         tempclm,tempplt,tempsnow,tempagg, &
     scmixg,scmixh
!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!         tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
       real refreclm,refimclm,refreclice,refimclice,refreagg, &
          refimagg,refregr,refimgr,refreha,refimha, &
          refreair,refimair,refreplt,refimplt
      

    real*8 aovrb
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!        (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!        (glbdy(7),deq),(glbdy(8),dmx), &
!        (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!        (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!        (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!

      refre=refreclice
      refim=refimclice
      scshp= 1    ! 1 denotes sphere/spheroid shape
      dmx=dsz
    
!      aovrb=1.0
      aovrb=.2
      deq=dmx/aovrb**(2.d0/3.d0)*min(1.0d0,aovrb)
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      if(dsz.gt. 1.0) cdtyp(1)=0
      canum(1)= 10
      calow(1)=  0
      caupp(1)= 10
      caavr(1)=  0
      cadev(1)= 5

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10
      calow(2)=  0
      caupp(2)=360
!
!

      return
      end subroutine 


             subroutine hydcldrp(dsz,epswtr,aovrb)


    implicit none
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     cloud droplet of characterictic size dsz=deq.
!-----------------------------------------------------------------------
!
      real dsz, glbdy(12), glcad(2,6), &
           tempenv,temprain,tempgraup,temphail,tempcld,&
           tempclm,tempplt,tempsnow,tempagg, &
           scmixg,scmixh
!           snowdens,aggdens,gradens,haildens,clmdens,pltdens
      complex epswtr, epsice
!      common/glhydbdy/glbdy,glcad
!      common/temperature/tempenv,temprain,tempgraup,temphail,tempcld,&
!           tempclm,tempplt,tempsnow,tempagg
!      common/heatbudget/scmixg,scmixh
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2), &
!           sctmp,scmix
      real*8 freq,aovrb
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!         (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!         (glbdy(7),deq),(glbdy(8),dmx), &
!         (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!         (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!         (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!

      
      call watereps(sctmp,epswtr,epsice)
      refre= real(sqrt(epswtr))
      refim=aimag(sqrt(epswtr))
    
      scshp= 1    !1 denotes sphere/spheroid shape
      deq=dsz
      aovrb=1.0
      dmx=deq*aovrb**(2.d0/3.d0)/min(1.0d0,aovrb)
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: one point Uniform distribution
!
      cdtyp(1)= 1
      if(dsz.gt.0.025) cdtyp(1)=0
      canum(1)= 1
      calow(1)=-1
      caupp(1)= 1

!-->  azimuthal angle: one point Uniform distribution
      cdtyp(2)= 1
      canum(2)= 1
      calow(2)= -1
      caupp(2)= 1


      return
      end subroutine

      subroutine hydrains(dsz,epswtr,aovrb)
!
!-----------------------------------------------------------------------
!     This routine defines particle information necessary to scattering
!     computation and canting angle distribution integration, for single
!     rain drop of characteristic size dsz=deq.
!-----------------------------------------------------------------------
!
    

    Implicit none

      complex epswtr, epsice

      integer, parameter ::mst=9
      real temp, scmixg, scmixh, glbdy(12), glcad(2,6),&
!    snowdens,aggdens,gradens,haildens,clmdens,pltdens, &
        tempenv,temprain,tempgraup,temphail,tempcld,&
        tempclm,tempplt,tempsnow,tempagg, dsz
      
!      common/glhydbdy/glbdy,glcad
      
!      common/heatbudget/scmixg,scmixh
!      common/density/snowdens,aggdens,gradens,haildens,clmdens,pltdens
!      real cdtyp(2),canum(2),calow(2),caupp(2),caavr(2),cadev(2), &
!         sctmp,scmix
      real*8 freq, aovrb
!
!-->  Global arrays glbdy, glcad are equivalent to these parameters:
!
!      equivalence (glbdy(1),freq),(glbdy(2),sctmp),(glbdy(3),scmix), &
!      (glbdy(4),refre),(glbdy(5),refim),(glbdy(6),scshp), &
!      (glbdy(7),deq),(glbdy(8),dmx), &
!     (glbdy(10),wcnxl),(glbdy(11),wcnzc),(glbdy(12),wcnlm)
!      equivalence (glcad(1,1),cdtyp(1)),(glcad(1,2),canum(1)), &
!       (glcad(1,3),calow(1)),(glcad(1,4),caupp(1)), &
!       (glcad(1,5),caavr(1)),(glcad(1,6),cadev(1))
!
!===> Define particle size, shape and electromagnetic properties.
!
      
!      call watereps(sctmp,epswtr,epsice)
      

      refre= real(sqrt(epswtr))
      refim= aimag(sqrt(epswtr))
      

      scshp= 1    ! 1 denotes sphere/spheroid shape
      deq=dsz
      

!
!     Set aspect ratio to one to speed up code for RAMS
!     processing
!
!      aovrb=1.0
!
!*********************************************************
!-->  Bear and Chung's aspect ratio formulae
!     Slight change in Beard/Chuang formula. Is from 
!     Andsager etaal   Bringi 
!
     if(deq.gt.1.) then        !vorher 0.3 mal zum testen raufgesetzt
       aovrb=1.0048+.0057*(deq/10.)-2.628*(deq/10.)**2 &
       +3.682*(deq/10.)**3-1.677*(deq/10.)**4
     else
       aovrb=1.0
     endif
!   
!*********************************************************
!     This is a Prupacher-Pitter formula : Equilibrium axis ratio
!mh   if (deq .gt. 0.49) then
!mh     aovrb = 1.03 - 0.062*deq
!mh   else
!mh     aovrb = 1.00
!mh   endif
 
!*********************************************************
!vorher gelaufen, jetzt 1. variation
!
!   Added Andsager fit for accounting for drop oscillations
!   Bringi 
      
!      if (deq.le.1..or.deq.gt.4.) then
!    Use Beard and Chuang equilibrium shapes
!        aovrb = 1.0048+0.0057*(deq/10.)-2.628*(deq/10.)**2+3.682*(deq/10.)**3-1.677*(deq/10.)**4
!      else
!     this Andsager and Beard fit for oscillating drops using
!     27 data points from Chandra et al,Beard,Kubesh and Beard
!     and 4 data points from Andsager et al.
        aovrb=1.012-0.1445*(deq/10.)-1.028*(deq/10.)**2
!      endif 
      

!**********************************************************
!  Use NR fit proposed by Keenan et al. (1999)
!mh   aovrb = 0.9939 + 0.00736*deq - 0.018485*deq*deq + 
!mh  &        0.001456*deq*deq*deq
!mh   if (aovrb .gt. 1.0) aovrb = 1.0
!**********************************************************
!
!===> Define orientation parameters of the body symmetry axis.
!
!-->  polar angle: Gaussian distribution
!
      cdtyp(1)= 2
      canum(1)= 10
      calow(1)=  0
      caupp(1)=  10
! caavr is the mean polar canting angle
      caavr(1)=  0
!***************************************************
!     changed std dev of canting angle to 10 deg 
!     for tropical rain. 
      cadev(1)= 5

!-->  azimuthal angle: Uniform distribution
      cdtyp(2)= 1
      canum(2)= 10 
      calow(2)= 10
      caupp(2)= 360

!  for single-drop use no orientation
!      if (glszd(2,1) .eq. 3) then
!!        cdtyp(1)= 1
!        canum(1)= 1
!        calow(1)= -1
!        caupp(1)= 1
!        caavr(1)= 0
!        cadev(1)= 0
!        cdtyp(2)= 1
!        canum(2)= 1
!        calow(2)= -1
!        caupp(2)= 1
!      endif


      return
      end subroutine

      end module src_hyd
