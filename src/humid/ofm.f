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
      subroutine ofm (kk,laps_p,laps_t,laps_q,laps_sfc_t,
     1     psfc, jday, lat, ZA,
     1                                Tbest, 
     1     sec_za,              !optran 90 variable 
     1     sfc_emis,            !optran 90 variable
     1     sfc_refl,            !optran 90 variable
     1     sec_solar            !optran 90 variable
     1     )


c   This routine interfaces GOES 8/10 satellite broadcast network data (and
c   local GVAR data) to the LAPS moisture analysis.  In 1999, this routine
c   was modified from an earlier version that used the University of
c   Wisconsin -- Madison's forward model to a newer model developed at
c   NESDIS.  OPTRAN (optical transmittance) forward model was developed by
c   Thomas Kleespies (NESDIS) and questions about this model should be
c   directed to him.  Forecast Systems Laboratory does not in any way
c   guarantee the validity of OPTRAN and distributes this software on an
c   as-is basis.  MOREOVER, FSL HAS PERMISSION TO DISTRIBUTE OPTRAN AS PART
c   OF LAPS TO FEDERAL AGENCIES.  NON-FEDERAL ENTITIES NEED TO INQUIRE WITH
c   NESDIS TO ESTABLISH THEIR RIGHTS AND OBLIGATIONS WITH REGARD TO OPTRAN.
c   
c   The version of OPTRAN with which this software is used, has been
c   modified by FSL to include both sounder and imager channels for a
c   particular satellite in one call to the routine.  Thus a user only need
c   to setup OPTRAN for a particular satellite.  After doing such, either
c   the imager or sounding instrument can be used with the software without
c   further recompilation.  
      
 

      USE module_sfc_structure

      Implicit None

      save
      
      Include '../include/Constants.inc'
      
      include '../include/trans.inc'
      
      integer Mchan
      parameter (Mchan=10)
      
      integer kk                ! kk the laps profile dim, 
      integer Nk                ! the size of the composite vectors
      integer start_level       !lowest level of climo used
      real laps_p(kk),laps_t(kk),laps_q(kk), lat, laps_sfc_t, psfc
      real ZA                   !zenith angle (degrees)
      real sec_za               !secant of ZA (for optran 90)
      real sfc_emis             !surface emissivity for optran 90
      real sfc_refl             !surface reflectivity for optran 90
      real sec_solar            !local secant of solar angle for optran 90
      integer jday, istatus
      logical first_call
      data first_call /.true./

c     optran 90 variables

      real, dimension (0:Nlevel,1)::  level_p,level_t, level_w, level_o
      real, dimension (Nlevel,1)  ::  layer_p,layer_t, layer_w,layer_o
      integer nk_90 ! optran 90 reduced nk counter
      
c     climo variables
      real c_p (40)
      real c_t (40)
      real c_m (40)
      
      real conversion
      real mdf
      
c     satellite number

      integer goes_number
      common /sat_id/ goes_number
      
      Real TotTrans(0:Nlevel,Nchan)
      
      Real O(Nlevel)
      
      Real Tbest(Nchan)         ! Brightness temps from estimated tau
      
      Real Compbright           ! brightness temp function
      
      Integer Ichan,Iatm,Iangle,Level
      
      Real Angle(5) 
      data angle /0.0,36.869898,48.189685,55.150095,60.0/
      
      integer Channels(Mchan)
      data channels     /10, 8, 7, 11, 16, 6, 12, 20, 21, 22 /
      
c     Integer Channels(Nchan) 
c     &       	/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 /
      
c     Real Vnu(Nchan)
      
c     goes sounder vnu
c     data vnu /680.,696.,711.,733.,748.,790.,832.,907.,1030.,
c     1           1345.,1425.,1535.,2188.,2210.,2248.,2420.,2513.,2671./
      
      integer Jatm,Ist,Ind,i
      
c     goes sounder vnu goes imager vnu
      
      real vnu(Mchan)
      
      data vnu / 1345., 907., 832., 1425., 2420., 790., 1535.,
     1     1482., 935., 835.4/


c     ----------------------do ozone first call ----

      if (first_call) then
         first_call = .false.
         
c     determine climo start level
         
         call climate_sm(lat, jday, c_p, c_t, c_m, istatus)
         
         start_level = 40       ! (high c_p)
         do while (  c_p(start_level) .ge. laps_p(kk) )
            start_level = start_level - 1
         enddo
         
c     fill pressure level
         
         nk = 0
         
         do i = 1,start_level
            nk = nk +1
            p(nk)   = c_p(i)
         enddo
         
         do i = kk,1,-1
            nk      = nk+1
            p(nk)   = laps_p(i)
         enddo
         
         write(6,*) 'Ozone levels used = ',nk
         
c     grab ozone
         
         call oh_zone (p,o,nk,1,istatus)
         
c     convert (g/g) to (g/kg)
         
         conversion    = 1.0e+3
         
c     ozone called only one time
         
      endif
      
c------------------------------------------------
      
c     00000000000000000000000000000
c     assemble vectors module 1
c000000000000000000000000000000
      
c     build vectors, count total levels (not including surface)
      
c     determine climo start level
      
      start_level = 40          ! (high c_p)
      do while (  c_p(start_level) .ge. laps_p(kk) )
         start_level = start_level - 1
      enddo

      nk = 0
      
      do i = 1,start_level
         nk = nk +1
         p(nk)   = c_p(i)
         t(nk,1) = c_t(i)
         q(nk,1) = c_m(i)       ! note this is g/kg mixing ratio
      enddo
      
      do i = kk,1, -1
         if(laps_p(i) .lt. psfc) then ! accept
            nk      = nk+1
            p(nk)   = laps_p(i)
            t(nk,1) = laps_t(i)
            q(nk,1) = laps_q(i) * conversion ! note this is g/kg mixing ratio
         endif
      enddo
      
c     now at sfc or cloud top
      
      nk = nk + 1
      p(nk) = psfc
      t(nk,1) = laps_sfc_t
      q(nk,1) = q(nk-1,1)       ! approximation for now
      
      if(nk.gt.Nlevel) then
         write(6,*) 'Array dimension error'
         write(6,*) 'Module ofm.f'
         write(6,*) 'Parameter Nlevel too small'
      endif
      
c0000000000000000000000
c     end assembling vectors
c00000000000000000000000

cc insert here new code for OPTRAN 90 INTERFACE  dB  2002
c upper level layer computations.

      do i = 1,nk
         level_p (i-1,1) = p(i)
         level_t (i-1,1) = t(i,1)
         level_w (i-1,1) = q(i,1)
         level_o (i-1,1) = o(i)
      enddo
      nk_90 = nk - 1            ! (start counter at zero for layers)

c     generate layers

      call make_optran_layers ( level_t, level_w, level_o, 
     1     layer_t, layer_w, layer_o, level_p,layer_p,
     1     nk_90, istatus)

c surface parameters

      


cc end insert new code
      
c11111111111111111111
c     begin using vectors
c11111111111111111111
      
      Iatm = 1 
      
      Call Precip_water_Profile(q(1,Iatm),P,0.0,Nk,Pw(1,Iatm))
      Call Precip_water_Profile(O        ,P,0.0,Nk,O3(1,Iatm))
      
      Iangle = 1
      Call Optrans(
     &     nk,
     &     P,
     &     T(1,Iatm),
     &     Q(1,Iatm),
     &     Pw(1,Iatm),
     &     O3(1,Iatm),
     &     ZA,
     &     Mchan,
     &     Channels,
     &     TotTrans)
      
c1111111111111111111
c     end using vectors
c1111111111111111111
      
c2222222222222222222
c     begin computing brightness
c2222222222222222222
      
      Do Ichan = 1 , Mchan
         
         TbEst(Ichan) = Compbright(Vnu(Ichan),T(1,Iatm),
     &        TotTrans(0,Ichan),laps_sfc_t,Nk)    
         
      EndDo
      
c2222222222222222222222
c     end btemp computation
c2222222222222222222222
      
      return
      
      End
