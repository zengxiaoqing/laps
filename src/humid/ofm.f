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
      subroutine ofm (kk,laps_p,laps_t,laps_q,laps_sfc_t,
     1     psfc, jday, lat, ZA,
     1                                Tbest)


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
      integer jday, istatus
      logical first_call
      data first_call /.true./
      
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
      
      Real Angle(5) /0.0,36.869898,48.189685,55.150095,60.0/
      
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
