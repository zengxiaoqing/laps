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
        Subroutine Normalize_Brightness(i4time,lat,lon,image,ni,nj
     1  ,sublat_d,sublon_d,range_m,l_national,iskip_bilin,r_missing_data
     1  ,lun,i_dir,Phase_angle_d,Specular_ref_angle_d,istatus)


C***Normalize a vis satellite image for solar angle.

C       S. Albers          Feb 94       Original version
C       J. Smart           Jan 96       Output Phase angle array
C       J. Smart           Mar 97       Added r_missing_data to argument list.
C       S. Albers          May 97       No longer using 0 as missing data
C                                       value. Added error trapping when
C                                       Emission angle < 0.
C       S. Albers          Jun 97       Error handling for 'iskip_bilin'

C Argument      I/O       Type                    Description
C --------      ---       ----    ------------------------------------------------
C i4time         I        I*4     i4time of image (1960 reference)
C lat            I        R*4  A  latitude array of image pixels (degrees)
C lon            I        R*4  A  longitude array of image pixels (degrees)
C image         I/O       R*4  A  Remapped image array.
C ni             I        I*4     i dimension of image
C nj             I        I*4     j dimension of image
C sublat_d       I        R*4     Latitude of satellite subpoint (degrees)
C sublon_d       I        R*4     Longitude...                   (degrees)
C range_m        I        R*4     Distance of spacecraft from center of earth (m)
C l_national     I        L       .true. for CONUS / .false. for LAPS
C iskip_bilin    I        I*4     >5 for CONUS / ~1 for LAPS
C                                 Controls subsampling of the grid for the more
C                                 time-consuming calculations. Larger values
C                                 give faster runtimes and less accuracy. 
C                                 Smaller values give slower runtimes and
C                                 greater accuracy, and also risks an error
C                                 message about dimensions. In general, the 
C                                 product of the satellite data resolution and
C                                 'iskip_bilin' should lie between about 10km 
C                                 and 30km.
C lun            I        I*4     Logical Unit # for logging output
C i_dir          I        I*4     Direction of normalization (-1,0,+1)
C phase_angle_d  O        R*4  A  Phase Angle (sparse array if iskip_bilin > 1)
C Specular_ref_angle_d O  R*4  A  Distance from specular reflection pt to sun
      include 'trigd.inc'
C istatus        O        I*4     Standard status return.

      include 'trigd.inc'
        integer         ni,nj
C***Parameter list variables
        real*4          lat(ni,nj),lon(ni,nj)
        Real*4          sublat_d,sublon_d,range_m
        Integer*4       i4time,istatus
        real*4          image(ni,nj)
        logical         l_national

C***Local variables
        integer*4 maxlut
        real*4 rpd,radius_earth_m
        Parameter       (maxlut=205,         ! largest array size allowed
     1          rpd=3.1415926536/180.,
     1          radius_earth_m = 6378137.)

        Real*8  TX,TY,TZ,RX,RY,RZ,
     1          SATX,SATY,SATZ,R8Phase_angle_r,R8Emission_angle_r,
     1          R8_ref_angle_r,
     1          stx_n,sty_n,stz_n,arg_dum,
     1          tx_n,ty_n,tz_n,
     1          refx,refy,refz

        Real*4  normfac,imgtmp,
     1          solar_factor(maxlut,maxlut),solar_alt_d,XFrac,YFrac,
     1          SF_UL,SF_UR,SF_LR,SF_LL,SF_U,SF_L,S_F,Weight,
     1                PF_UL,PF_UR,PF_LR,PF_LL,PF_U,PF_L,P_F,RBril,RBrih,
     1          Phase_factor(maxlut,maxlut),Phase_angle_d(ni,nj),
     1          sat_radius,Emission_angle_d,Specular_ref_angle_d(ni,nj)       

        Integer nilut,njlut,I,J,img_i(maxlut),img_j(maxlut),
     1            ilut,jlut,ISpace,JSpace, ni2, nj2


        write(lun,*)' Begin normalization routine (V970519)'
        write(lun,*)' Sub Lat/Lon/Rng = ',sublat_d,sublon_d,range_m

        iwrite = 0

C***Where's the sun?
!       Use cartesian coordinates with the x-axis through the prime
!       meridian and distances in AU.
        rlat = 0.
        rlon = 0.
        call solar_position(rlat,rlon,i4time,solar_alt_d,solar_dec_d,hr_
     1angle_d)
        solar_range = 1.
        solar_sublon_d = -hr_angle_d
        RX = cosd(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        RY = sind(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        RZ = sind(solar_dec_d)                        * solar_range

!   Satellite Location (in AU)
        au_m = 149000000.
        sat_radius = range_m / au_m

        SATX = cosd(sublon_d) * cosd(sublat_d) * sat_radius
        SATY = sind(sublon_d) * cosd(sublat_d) * sat_radius
        SATZ = sind(sublat_d)                  * sat_radius

        write(lun,*)'    I    J    ALT    EMIS   PHA    PF   VIS  SPEC'       

C***Fill the solar brightness and phase angle arrays
        normfac=sind(58.)       ! Normalized sun angle

        nilut = (ni-2) / iskip_bilin + 2
        njlut = (nj-2) / iskip_bilin + 2

        if(nilut .gt. maxlut .or. njlut .gt. maxlut)then
            write(lun,*)'Insufficient dimension for maxlut'
            write(lun,*)'maxlut = ',maxlut
            write(lun,*)'nilut = ',nilut
            write(lun,*)'njlut = ',njlut

            rskip_new = max(float(nilut)/float(maxlut)
     1                     ,float(njlut)/float(maxlut))
            iskip_new = int(rskip_new) + 1

            write(lun,*)'Input parameter iskip_bilin currently equals'   
     1                 ,iskip_bilin
            write(lun,*)'Try increasing iskip_bilin to approx '
     1                 ,iskip_new      
            write(lun,*)'Secondary alternative is to increase maxlut'
     1                 ,' declared as a parameter in normalize.f'
            write(lun,*)
     1      'Setting istatus to zero, returning without normalizing'

            istatus = 0
            return
        else
            write(lun,*)' Dimension for maxlut',maxlut,nilut,njlut
        endif
        nj2 = nj
        ni2 = ni
        Do jlut=1,njlut
         j = ((jlut-1) * iskip_bilin) + 1
         
         j = min(j,nj2)
         img_j(jlut) = j

         Do ilut=1,nilut

          i = ((ilut-1) * iskip_bilin) + 1
          i = min(i,ni2)
          img_i(ilut) = i

C   Compute equatorial coordinates of point on sfc of earth
          TX = cosd(lon(i,j)) * cosd(lat(i,j)) * radius_earth_m / au_m
          TY = sind(lon(i,j)) * cosd(lat(i,j)) * radius_earth_m / au_m
          TZ = sind(lat(i,j))                  * radius_earth_m / au_m

          call solar_position(lat(i,j),lon(i,j),i4time,solar_alt_d
     1                                  ,solar_dec_d,hr_angle_d)

C   Reduce and limit correction at terminator
          If(solar_alt_d.lt.8.)then
           solar_alt_d=max(8.-(8.-solar_alt_d)*.5,3.8)
          EndIf

          solar_factor(ilut,jlut)=log(sind(solar_alt_d)/normfac)

C   Compute Emission Angle (Emission_angle_d = satellite angular altitude)
          Call AngleVectors(SATX-TX,SATY-TY,SATZ-TZ,TX,TY,TZ
     1                            ,R8Emission_angle_r)
          Emission_angle_d = 90. - R8Emission_angle_r / rpd

          if(Emission_angle_d .lt. 0.)then
              Emission_angle_d = 0.
              istatus = 0

              if(iwrite .le. 1)then
                If(image(i,j) .ne. r_missing_data)then ! valid image data
                  write(6,*)
     1            ' Warning, Emission_angle_d =', Emission_angle_d
                  write(6,*)
     1            ' You are normalizing beyond the earths limb'     
                  write(6,*)
     1            ' Check your satellite subpoint and lat/lons'
                  write(6,*)'i,j,lat(i,j),lon(i,j),sublat_d,sublon_d'        
                  write(6,*)i,j,lat(i,j),lon(i,j),sublat_d,sublon_d
                  iwrite = iwrite + 1
                endif ! valid image pixel intensity
              endif ! iwrite .le. 1

          endif ! emission angle < 0 (looking beyond the limb)

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


          if(i .lt. 40)then ! for testing

          Phase_factor(ilut,jlut) = cosd(phase_angle_d(i,j))**6
!    1                    * sind(emission_angle_d)**1.0
!    1                    * sind(solar_alt_d)**1.0
     1                    * sind(emission_angle_d)**0.5
     1                    * sind(solar_alt_d)**0.25

          else
          Phase_factor(ilut,jlut) = cosd(phase_angle_d(i,j))**6
!    1                    * sind(emission_angle_d)**1.0
!    1                    * sind(solar_alt_d)**1.0
     1                    * sind(emission_angle_d)**0.5
     1                    * sind(solar_alt_d)**0.125


          endif


!         Phase_factor(ilut,jlut) = 0.    ! Disable Phase Factor Correction

!         Perhaps the Phase factor should be modified if the solar altitude
!         is low (< 23 deg). This would probably be designed to darken
!         the land areas at low phase angle and low solar altitude.


          if(ilut .eq. ilut/10*10 .and. jlut .eq. jlut/10*10)then

              imgtmp=image(i,j)
              if(imgtmp.eq.r_missing_data)imgtmp=-99.00
              write(lun,50)ilut,jlut,solar_alt_d
     1                  ,emission_angle_d,phase_angle_d(i,j)
     1                  ,phase_factor(ilut,jlut),imgtmp
     1                  ,specular_ref_angle_d(i,j)
     1                  ,lat(i,j),lon(i,j)
 50           format(1x,2i5,f7.2,f7.2,f7.2,f6.2,2f6.1,2f7.2,' __')      
          endif

         EndDo
        EndDo

        write(lun,*)' Normalization Lookup Tables complete'

C***Apply the solar brightness normalization to the image
        jlut=1
        JSpace=img_j(jlut+1)-img_j(jlut)

        Do j=1,nj

         If(j .gt. img_j(jlut+1))then
          jlut=jlut+1
          JSpace=img_j(jlut+1)-img_j(jlut)
         EndIf

         YFrac=Float(j-img_j(jlut))/JSpace

         ilut=1
         ISpace=img_i(ilut+1)-img_i(ilut)

         Do i=1,ni

          If(i .gt. img_i(ilut+1))then
           ilut=ilut+1
           ISpace=img_i(ilut+1)-img_i(ilut)
          EndIf

          XFrac=Float(i-img_i(ilut))/ISpace

C   Bilinearly interpolate the normalization factor from the surrounding points.
          SF_UL=solar_factor(ilut  ,jlut  )
          SF_UR=solar_factor(ilut+1,jlut  )
          SF_LR=solar_factor(ilut+1,jlut+1)
          SF_LL=solar_factor(ilut  ,jlut+1)

          SF_U=SF_UL+(SF_UR-SF_UL)*XFrac
          SF_L=SF_LL+(SF_LR-SF_LL)*XFrac
          S_F=SF_U+(SF_L-SF_U)*YFrac

C   Bilinearly interpolate the phase factor from the surrounding points.
          PF_UL=Phase_Factor(ilut  ,jlut  )
          PF_UR=Phase_Factor(ilut+1,jlut  )
          PF_LR=Phase_Factor(ilut+1,jlut+1)
          PF_LL=Phase_Factor(ilut  ,jlut+1)

          PF_U=PF_UL+(PF_UR-PF_UL)*XFrac
          PF_L=PF_LL+(PF_LR-PF_LL)*XFrac
          P_F=PF_U+(PF_L-PF_U)*YFrac

C   Greater (lesser) abs. values of RBriH and RBriL will brighten (darken) the
C   high and low ends, respectively.  RBriL is modified to make darker land more
C   uniform in brightness.
          RBriH=-60.*S_F   -   20.*P_F

          if(l_national)then
            Weight=18.
          Else ! Local type scales
            If(S_F.gt.-1.0)then
             Weight=30.
            Else
             Weight=Max(30.-12.*(-1.0-S_F),20.)
            EndIf
          EndIf

          RBriL=-Weight*S_F  -   20.*P_F

          If(image(i,j) .ne. r_missing_data)then
            if(i_dir .eq. +1)then
              Call Stretch(68.-RBriL,220.-RBriH,68.,220.,image(i,j))
            elseif(i_dir .eq. -1)then
              Call Stretch(68.,220.,68.-RBriL,220.-RBriH,image(i,j))
            endif
          endif

         EndDo
        EndDo

        istatus = 1
        write(lun,*)' Normalization complete'

        Return

        End
C-------------------------------------------------------------------------------
        Subroutine Stretch(IL,IH,JL,JH,rArg)

        Implicit        None

        Real*4          A,B,IL,IH,JL,JH,rarg

        a = (jh - jl) / (ih - il)
        b =  jl - il * a

        rarg = a * rarg + b
        rarg = max(min(255.,rarg),0.)

        return
        end


        subroutine deviate_ray(l1,x5,y5,z5,x3,y3,z3,x8,y8,z8)

	implicit real*8(a-z)

        c9 = x3*x5 + y3*y5 + z3*z5
        c8 = sqrt(1.d0 - l1**2 * (1.d0 - c9**2))
        m1 = c8 - l1*c9

        x8 = l1*x5 + m1*x3                    ! x dir cosine of outgoing ray
        y8 = l1*y5 + m1*y3                    ! y dir cosine of outgoing ray
        z8 = l1*z5 + m1*z3                    ! z dir cosine of outgoing ray

	return
	end



