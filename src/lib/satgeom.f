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
        Subroutine satgeom(i4time,lat,lon,ni,nj
     1  ,sublat_d_a,sublon_d_a,range_m,r_missing_data,Phase_angle_d
     1  ,Specular_ref_angle_d,Emission_angle_d,azimuth_d,istatus)


        include 'trigd.inc'

        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.

C istatus        O        I*4     Standard status return.
        integer         ni,nj
C***Parameter list variables
        real          lat(ni,nj),lon(ni,nj)
        real          sublat_d_a(ni,nj),sublon_d_a(ni,nj)
        Real          sublat_d,sublon_d,range_m ! sat from Earth center
        Integer       i4time,istatus

C***Local variables
        integer maxlut
        real rpd
        Parameter   (maxlut=2500,         ! largest array size allowed
     1          rpd=3.1415926536/180.)

        Real*8  TX,TY,TZ,RX,RY,RZ,
     1          SATX,SATY,SATZ,R8Phase_angle_r,
     1          R8_ref_angle_r,
     1          stx_n,sty_n,stz_n,arg_dum,
     1          tx_n,ty_n,tz_n,
     1          refx,refy,refz

! Note that solar_factor and phase_factor originally were declared to maxlut
        Real  normfac,imgtmp,                                                                               
     1          solar_alt_d(ni,nj),
     1          XFrac,YFrac,
     1          SF_UL,SF_UR,SF_LR,SF_LL,SF_U,SF_L,S_F,Weight,
     1                PF_UL,PF_UR,PF_LR,PF_LL,PF_U,PF_L,P_F,RBril,RBrih,
     1          Phase_factor(ni,nj),Phase_angle_d(ni,nj),
     1          sat_radius,Emission_angle_d(ni,nj),
     1          azimuth_d(ni,nj),
     1          Specular_ref_angle_d(ni,nj)       

        Real RBril_a(ni,nj),RBrih_a(ni,nj)

        Integer nilut,njlut,I,J,img_i(maxlut),img_j(maxlut),
     1            ilut,jlut,ISpace,JSpace, ni2, nj2

        character*9 a9time

        call make_fnam_lp (i4time,a9time,istatus)
        if(istatus .ne. 1)return

        lun = 6
        write(lun,*)' Begin satgeom at ',a9time       

        call zero(phase_angle_d,ni,nj)
        call zero(emission_angle_d,ni,nj)
        call zero(specular_ref_angle_d,ni,nj)

        iwrite = 0

C***Where's the sun?
!       Use cartesian coordinates with the x-axis through the prime
!       meridian and distances in AU (geocentric equatorial).
        rlat = 0.
        rlon = 0.
        call solar_position(rlat,rlon,i4time,solar_alt_deg       
     1                     ,solar_dec_d,hr_angle_d)
        solar_range = 1.
        solar_sublon_d = -hr_angle_d
        RX = cosd(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        RY = sind(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        RZ = sind(solar_dec_d)                        * solar_range

!   Satellite Location (in AU - geocentric equatorial)
        au_m = 149000000.
        sat_radius = range_m / au_m

        write(lun,*)'    I    J    ALT    EMIS   PHA    PF   VIS  SPEC'       

C***Fill the solar brightness and phase angle arrays
        normfac=sind(58.)       ! Normalized sun angle

        Do j = 1,nj
         Do i = 1,ni

          sublat_d = sublat_d_a(i,j)
          sublon_d = sublon_d_a(i,j)

          SATX = cosd(sublon_d) * cosd(sublat_d) * sat_radius
          SATY = sind(sublon_d) * cosd(sublat_d) * sat_radius
          SATZ = sind(sublat_d)                  * sat_radius

          intvl = max(ni/21,1)
          if(j .eq. nj/2 .AND. i .eq. (i/intvl)*intvl)then
              idebug = 1
          else
              idebug = 0
          endif

C   Compute Emission Angle (Emission_angle_d = satellite angular altitude)
 
          call sat_angular_alt(sat_radius,lat(i,j),lon(i,j)
     .,SATX,SATY,SATZ,TX,TY,TZ,Emission_angle_d(i,j),istatus)

          if(Emission_angle_d(i,j) .lt. 0.)then
              istatus = 0
              Emission_angle_d(i,j) = 0.0
          endif ! emission angle < 0 (looking beyond the limb)

!         Topocentric equatorial vector of satellite relative to gridpoint
          DX=SATX-TX
          DY=SATY-TY
          DZ=SATZ-TZ

          if(idebug .eq. 1)then
              write(6,*)
              write(6,*)'/lon/sub/DX/DY/DZ',lon(i,j),sublon_d,DX,DY,DZ       
          endif

!         Rotate this vector around Z axis to get local cartesian coordinates
          call rotate_z(DX,DY,DZ,-lon(i,j))

!         Convert cartesian coordinates to dec and ha
          call xyz_to_polar_d(DX,DY,DZ,dec,angle,r)

          if(angle .gt. 180.)then
              angle = angle - 360.
          endif

          ha = -angle

          if(idebug .eq. 1)then
              write(6,*)'rotated DX/DY/DZ,dec,ha',DX,DY,DZ,dec,ha
          endif

!         Convert dec and ha to alt/az
          call equ_to_altaz_d(dec,ha,lat(i,j),alt,azimuth_d(i,j))

          if(idebug .eq. 1)then
              write(6,*)'alt/az',alt,azimuth_d(i,j)
          endif

          goto500

          call solar_position(lat(i,j),lon(i,j),i4time
     1                                  ,solar_alt_d(i,j)
     1                                  ,solar_dec_d,hr_angle_d)

C   Compute Phase Angle
          Call AngleVectors(TX-SATX,TY-SATY,TZ-SATZ,RX,RY,RZ
     1                     ,R8Phase_angle_r)
          Phase_angle_d(i,j) = 180. - R8Phase_angle_r / rpd

C   Compute Specular Reflection Angle
          stx_n = tx-satx
          sty_n = ty-saty
          stz_n = tz-satz
          call normalize(stx_n,sty_n,stz_n,arg_dum)

          tx_n = -tx
          ty_n = -ty
          tz_n = -tz
          call normalize(tx_n,ty_n,tz_n,arg_dum)

          call deviate_ray(-1d0,stx_n,sty_n,stz_n,tx_n,ty_n,tz_n
     1                         ,refx,refy,refz)
          Call AngleVectors(refx,refy,refz,RX,RY,RZ,R8_ref_angle_r)
          Specular_ref_angle_d(i,j) = 180. - R8_ref_angle_r / rpd

500       continue

         EndDo
        EndDo

c-----------------------------------------------------------------------
c=======================================================================

        return
        end

