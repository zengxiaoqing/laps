

       subroutine vi(magn_r8,elong_r8,altdif_r8,al1_r8,al2_r8,elgms_r8,elgmc_r8 &
                    ,al_t_r8,als_ct_r8,alm_best_r8,od_atm_a &
                    ,maglimd_r8,maglimt_r8,maglimn_r8,maglim_r8 &
                    ,vis_r8,c_observe)

       include 'trigd.inc'

       IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)

!      vis of 6.0 means barely visible to the naked eye, lower amounts are more visible
!      vis of 8.0 means barely visible to binoculars, lower amounts are more visible

       character*1 c_observe

       real*8 magn_r8, maglimd_r8, maglimt_r8, maglimn_r8, maglim_r8

       real   magn,elong,altdif,al1,al2,elgmc,elgms
       real   airmass,extinction,zenext,patm
       real   maglim,maglimd,maglimt,maglimm,maglimn,maglimp
       real   maglimm1,mag_eff
       real   al_t,als_ct,alm_best,sb1,sb2,rmaglim_to_b,b_to_maglim,maglimdark,sbdark
       real   al_t_eff,od_atm_a
       real   maglimdark2,sb_corr,diff_mag

       parameter (nmag_n=9)  ; real maglim_n(nmag_n); real alts_n(nmag_n)
       parameter (nmag_m=8)  ; real maglim_m(nmag_m); real alts_m(nmag_m)
       parameter (nmag_t=6)  ; real maglim_t(nmag_t); real alts_t(nmag_t)
       parameter (nmag_d=19) ; real maglim_d(nmag_d); real alts_d(nmag_d)

       magn = magn_r8
       elong = elong_r8
       altdif = altdif_r8
       al1 = al1_r8
       al2 = al2_r8
       elgms = elgms_r8
       elgmc = elgmc_r8
       al_t = al_t_r8
       als_ct = als_ct_r8
       alm_best = alm_best_r8

       maglimdark = 6.0
       sbdark = rmaglim_to_b(maglimdark)

       patm = 0.85
       call calc_extinction(90.,patm,airmass,zenext)

!      Nighttime (based on al1) http://www.asterism.org/tutorials/tut28-1.htm

!      Moonlight (based on lunar elongations for a full moon)
       i = 0
       i = i+1; alts_m(i) = 180. ;  maglim_m(i) = +4.8
       i = i+1; alts_m(i) = 90.  ;  maglim_m(i) = +5.3
       i = i+1; alts_m(i) = 45.  ;  maglim_m(i) = +4.4
       i = i+1; alts_m(i) = 10.  ;  maglim_m(i) = +3.8
       i = i+1; alts_m(i) =  5.  ;  maglim_m(i) = +3.3
       i = i+1; alts_m(i) =  2.  ;  maglim_m(i) = +2.0
       i = i+1; alts_m(i) =  1.  ;  maglim_m(i) = +1.0
       i = i+1; alts_m(i) =  0.25;  maglim_m(i) = -3.0

!      Twilight (based on altdif)
       i = 0
!      i = i+1; alts_t(i) = 28.   ;  maglim_t(i) =  5.0  ! consider extinction (alt = 12, value = 1.3)
!      i = i+1; alts_t(i) = 26.   ;  maglim_t(i) =  4.7  ! consider extinction (alt = 10, value = 1.6)
!      i = i+1; alts_t(i) = 23.   ;  maglim_t(i) =  3.9  ! consider extinction (alt =  7, value = 2.2) + skyglow 0.62
       i = i+1; alts_t(i) = 23.   ;  maglim_t(i) =  3.75 ! consider extinction (alt =  7, value = 2.2) + skyglow 0.94
!      i = i+1; alts_t(i) = 23.   ;  maglim_t(i) =  3.6  ! consider extinction (alt =  7, value = 2.2) + skyglow 1.29
       i = i+1; alts_t(i) = 16.   ;  maglim_t(i) =  2.4
       i = i+1; alts_t(i) = 12.   ;  maglim_t(i) =  0.0
       i = i+1; alts_t(i) =  8.   ;  maglim_t(i) = -3.0
       i = i+1; alts_t(i) =  3.   ;  maglim_t(i) = -6.0  ! consider extinction
       i = i+1; alts_t(i) =  0.001;  maglim_t(i) = -12.0 ! consider extinction
       if((elong) .lt. 0.5)then                          ! solar contrast effect
           if((elong) .ge. 0.25)then
               maglim_t(i) = maglim_t(i) - (20.*(0.5-elong))
           else
               maglim_t(i) = -17.0
           endif
       endif

!      Daylight (based on elong)
       i = 0
       i = i+1; alts_d(i) =180.  ;  maglim_d(i) = -2.75
       i = i+1; alts_d(i) =170.  ;  maglim_d(i) = -2.74
       i = i+1; alts_d(i) =160.  ;  maglim_d(i) = -2.68
       i = i+1; alts_d(i) =150.  ;  maglim_d(i) = -2.61
       i = i+1; alts_d(i) =140.  ;  maglim_d(i) = -2.50
       i = i+1; alts_d(i) =130.  ;  maglim_d(i) = -2.38
       i = i+1; alts_d(i) =120.  ;  maglim_d(i) = -2.24
       i = i+1; alts_d(i) =110.  ;  maglim_d(i) = -2.12
       i = i+1; alts_d(i) =100.  ;  maglim_d(i) = -2.03
       i = i+1; alts_d(i) = 90.  ;  maglim_d(i) = -2.0
       i = i+1; alts_d(i) = 80.  ;  maglim_d(i) = -2.1
       i = i+1; alts_d(i) = 30.  ;  maglim_d(i) = -3.0
       i = i+1; alts_d(i) = 20.  ;  maglim_d(i) = -3.3
       i = i+1; alts_d(i) = 12.  ;  maglim_d(i) = -3.7
       i = i+1; alts_d(i) =  9.  ;  maglim_d(i) = -4.0
       i = i+1; alts_d(i) =  2.  ;  maglim_d(i) = -6.0
       i = i+1; alts_d(i) =  1.  ;  maglim_d(i) = -8.0
       i = i+1; alts_d(i) =  0.5 ;  maglim_d(i) = -10.0
       i = i+1; alts_d(i) =  0.25;  maglim_d(i) = -17.0

       c_observe = ' '

       maglimn = -99.9
       maglimm = -99.9
       maglimp = -99.9
       maglimt = -99.9
       maglimd = -99.9
       maglim  = 99.9
       vis     = 99.9

       j = 0 ! debug flag

!      Nighttime limiting magnitude
       if(al1 .gt. 0.0)then
           call correct_skyglow(sbdark,al1,patm,sb2,diff_mag)
           maglimn = b_to_maglim(sb2)
           call calc_extinction(al1,patm,airmass,extinction)
           maglimn = maglimn + (zenext - extinction)
!          write(13,*)' al1/extinction/diff_mag/maglimn ',al1,zenext-extinction,diff_mag,maglimn
       endif

!      Moonlight limiting magnitude (when moon is up)
       if(al2 .lt. al1)then                        
           call interp_1d(alts_m,maglim_m,nmag_m,elgmc,maglimm1,-99.9,j)
           phase_angle_deg = 180. - elgms

!          Phase angle correction for moon's brightness
           mode = 1
           call phase_func_moon(phase_angle_deg,mode,area_rel,area_es,phase_corr)
           if(.true.)then                             
               sb1 = rmaglim_to_b(maglimm1) - sbdark
               sb2 = (sb1 / (10.**(phase_corr/2.5))) + sbdark
!              write(13,11)phase_corr,maglimm1,sb1,sb2,sbdark,b_to_maglim(sb2)
11             format(3x,' phase_corr: maglimm1/sb1/sb2/sbdark/b_to_maglim ',f8.2,3x,5f10.2)
               maglimm1 = b_to_maglim(sb2)
           endif

!          Sky brightness correction based on moon's altitude
           sb_corr = 2.0 * (1.0 - (sind(alm_best)**0.5))
           if(.true.)then                                      
               sb1 = rmaglim_to_b(maglimm1) - sbdark
               sb2 = (sb1 / (10.**(sb_corr/2.5))) + sbdark
!              write(13,12)sb_corr,maglimm1,sb1,sb2,sbdark,b_to_maglim(sb2)
12             format(3x,' sb_corr:    maglimm1/sb1/sb2/sbdark/b_to_maglim ',f8.2,3x,5f10.2)
               maglimm1 = b_to_maglim(sb2)
           endif

!          Combined sky brightness correction based on comet's altitude
           if(.true.)then                                      
               sb1 = rmaglim_to_b(maglimm1)                        
               call correct_skyglow(sb1,al1,patm,sb2,diff_mag)
!              write(13,13)maglimm1,sb1,sb2,b_to_maglim(sb2)
13             format(3x,' sb_corr:    maglimm1/sb1/sb2/b_to_maglim ',f8.2,3x,4f10.2)
               maglimm1 = b_to_maglim(sb2)
           endif

           maglimm = maglimm1 + (zenext - extinction) ! add extinction component
!          write(13,21)elgmc,phase_corr,alm_best,sb_corr,maglimm1,maglimn,maglimm
21         format(20x,' Moon is up: elgmc/phase_corr/alm_best/sb_corr/maglimm1/maglimn/maglimm = ',7f9.2)
       endif

!      Consider extinction/twilight based magnitude at moon rise/set                        
       if(al2 .lt. al1 .AND. al2 .gt. 0.)then                      
           if(al2 .ge. 7.0)then ! extinction at night
               call correct_skyglow(sbdark,al2,patm,sb2,diff_mag)
               maglimdark2 = b_to_maglim(sb2)
               call calc_extinction(al2,patm,airmass,extinction)
               maglimp = maglimdark2 + (zenext - extinction) ! add extinction component
           else ! twilight
               call interp_1d(alts_t,maglim_t,nmag_t,al2+16.,maglimp,-99.9,j)
           endif

!          write(13,22)maglimm,maglimn,maglimp
22         format(20x,' Moon rise/set: maglimm/maglimn/maglimp = ',68x,3f9.2)
       endif

!      Twilight limiting magnitude
       if(al1 .gt. 3.0)then
          call interp_1d(alts_t,maglim_t,nmag_t,al1+16.,maglimt,-99.9,j)
       else
          call interp_1d(alts_t,maglim_t,nmag_t,altdif,maglimt,-99.9,j)
       endif

!      Daylight limiting magnitude
       call interp_1d(alts_d,maglim_d,nmag_d,elong, maglimd,-99.9,j)
       if(al_t .gt. 0. .AND. maglimd .gt. -99.9)then ! correct for extinction
           call calc_extinction(al_t,patm,airmass,extinction)
           maglimd = maglimd + (zenext - extinction)

!          Correct for Sky Brightness change for solar altitude
           sb_corr = 2.0 * (1.0 - (sind(als_ct)**0.5))
           maglimd = maglimd + sb_corr

!          Correct for Sky Brightness change with comet altitude
           al_t_eff = al_t * min(max(1.0-(od_atm_a-0.05)*1.0, 0.0),1.0)
           call correct_skyglow(1.,al_t_eff,patm,sb2,sb_corr)
           maglimd = maglimd + sb_corr

!          write(13,31)elong,al_t,als_ct,sb_corr,zenext - extinction,maglimd
31         format(' Daylight - elong,al_t,als_ct,sb_corr,ext,maglimd: ',6f8.2)
       else
!          write(13,32)elong,al_t,als_ct,sb_corr,maglimd
32         format(' Skip Daylight - elong,al_t,als_ct,sb_corr,maglimd: ',5f8.2)
       endif

!      Apply for nighttime
       if(maglimn .gt. -99.9)then
!          call correct_skyglow(sbdark,al1,patm,sb2,diff_mag)
!          maglimn = b_to_maglim(sb2)
           maglim = maglimn
           vis = maglimdark + (magn - maglim)  
           c_observe = 'N'
!          write(13,*)' Correct nighttime mag for skyglow al1/maglimn = ',al1,maglimn
       endif

!      Apply for moonlight (when moon is up)
       if(c_observe .eq. 'N' .AND. maglimm .gt. -99.9 &
                             .AND. maglimm .lt. maglimn)then
           maglim = maglimm
           c_observe = 'M'
           vis = maglimdark + (magn - maglim)  
       endif

!      Apply for moonlight (extinction/twilight at moon rise/set time)  
       if(                         maglimp .gt. -99.9 &
                             .AND. maglimp .gt. maglim)then
           maglim = maglimp
           c_observe = 'O' ! Better at moon rise/set       
           vis = maglimdark + (magn - maglim)  
       endif

!      Apply for twilight
       if(maglimt .gt. -99.9)then
!          if(maglimt .lt. maglimm .OR. maglimm .eq. -99.9)then
               maglim = maglimt
               vis = maglimdark + (magn - maglim)  
               c_observe = 'T'
!          endif
       endif

!      Apply for daylight
       if(al_t .gt. 0. .AND. maglimd .gt. maglimt .AND. maglimd .gt. maglimn)then
!          write(13,*)' passed daylight test ',al_t,maglimd,maglimt
           vis_tmp = maglimdark + (magn - maglimd)  
           if(vis_tmp .le. 12.5)then
               maglim = maglimd
               vis = vis_tmp                    
               c_observe = 'D'
           endif
       endif

!      write(6,*)' magn,elong,altdif,maglim,vis,c_observe', &
!                  magn,elong,altdif,maglim,vis,c_observe

       maglimd_r8 = maglimd
       maglimt_r8 = maglimt
       maglimn_r8 = maglimn
       maglim_r8 = maglim
       vis_r8 = vis

       return
       end

!      Note Bortle relationship for solar elongation m = -2.5 + (2 * log(elong)) for 80mm scope
!                                           possibly m = -5.9 + (2 * log(elong)) for naked eye

!      Solar/Lunar altitude correction term
!      2 mags at zenith angle = 90.
!      0.5 mag at zenith angle = 30.
!      input solar/lunar altitude

       subroutine calc_extinction_simple(alt,airmass,extinction)

       airmass = 1. / sind(alt+1.4325)
       extinction = 0.28 * airmass

       return
       end

       subroutine calc_extinction(alt,patm,airmass,extinction)

       z = 90. - alt
       airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))
       extinction = 0.28 * airmass * patm

!      write(13,*)' calc_extinction: alt/airmass/extinction ',alt,airmass,extinction

       return
       end

       subroutine correct_skyglow(skyglow_in,alt,patm,skyglow_out,diff_mag)

!      Calculate change in skyglow from zenithal value for an idealized case

       mode = 1

       ratmax = 3.3
       if(alt .gt. 0.)then
           alt2 = alt + 10. * (1. - alt/90.)
           ratio = (1.0/sind(alt2))**0.68
           if(ratio .gt. ratmax)ratio = ratmax
       else
           ratio = ratmax
       endif

!      REM  Dark night sky brightness
!      BN=BO(I)*(1+.3*COS(6.283*(Y-1992)/11))    (Solar Cycle Term)
!      BN=BN*(.4+.6/SQR(1.0-.96*((SIN(ZZ))^2)))  (Zenith Angle Term)
!      BN=BN*(10^(-.4*K(I)*X))

!      ZZ = 90. - alt

!      BO = 1.0E-13
!      BN = skyglow_in
!      BN = BN*(.4+.6/SQRT(1.0-.96*((SIND(ZZ))**2)))
!      BN = BN*(10^(-.4*K(I)*X))

       if(mode .eq. 2)then
           elong = 90.
           call get_skyglow(0.,90.,90.,elong,patm,skyglowz)
           call get_skyglow(0.,90.,alt,elong,patm,skyglowo)
           ratio = skyglowo / skyglowz
       endif

       skyglow_out = skyglow_in * ratio

       diff_mag = -log10(ratio) * 2.5

       return
       end

       subroutine get_skyglow(rmag,altsource,altobj,elong,patm,skyglow)

!      REM  Daylight brightness
!      C4=10.0^(-.4*K(I)*XS)                        (Extinction and Airmass) 
!      FS=6.2E+07*(RS^-2)+(10^(6.15-RS/40))         (Elongation)
!      FS=FS+(10^5.36)*(1.06+((COS(RS*RD))^2))      (Elongation)
!      BD=10^(-.4*(MS(I)-MO(I)+43.27))              (Magnitude-Const)
!      BD=BD*(1-10^(-.4*K(I)*X))
!      BD=BD*(FS*C4+440000.0*(1-C4))

       real, parameter :: MO = -11.5

       call calc_extinction(altsource,patm,airmass,totexts)
       C4 = 10.0**(-.4*totexts)
       FS=6.2E+07*(elong**-2.0)+(10**(6.15-elong/40.))
       FS=FS+(10**5.36)*(1.06+((COSD(elong))**2))    ! Rayleigh Scattering
       BD=10.**(-.4*(rmag-MO+43.27))

       call calc_extinction(altobj   ,patm,airmass,totexto)
       BD=BD*(1.-10.**(-.4*totexto))
       BD=BD*(FS*C4+440000.0*(1.-C4))

       skyglow = BD/1.11E-15 ! convert to nanoLamberts

       return
       end
